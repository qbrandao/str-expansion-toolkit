#!/usr/bin/env python3
"""Per-sample subsampling fractions from measured depths.

Reads the collected output of 15_measure_depth.sbatch and emits one fraction
per sample to bring the whole panel to a single target depth. Nothing is
predicted here: every depth is measured, which is why this replaces the earlier
attempt to infer depth from file size.

The target defaults to the shallowest sample in the panel, so that no sample
has to be used above the depth at which every sample can be compared. Pass
--target to match the patient cohort depth instead, which is the usual choice.

Usage
  python3 make_downsample_plan.py measured_depth.tsv --target 20
  python3 make_downsample_plan.py measured_depth.tsv --target-from-min
"""
import argparse, csv, statistics, sys

NEEDED = ("label", "depth_x", "read_n50")
SPAN_SIZES = (1000, 2000, 5000, 10000, 20000, 50000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("measured")
    ap.add_argument("--target", type=float, default=None,
                    help="target depth in x, normally the patient cohort depth")
    ap.add_argument("--target-from-min", action="store_true",
                    help="use the shallowest sample in the panel as the target")
    ap.add_argument("--target-spanning", type=float, default=None,
                    help="target the expected number of reads spanning a locus of "
                         "--at-locus-bp instead of the raw depth. This is the better "
                         "target: two samples at the same raw depth can differ several "
                         "fold in spanning reads, because spanning depth depends on read "
                         "length as well as on depth.")
    ap.add_argument("--at-locus-bp", type=int, default=1000,
                    choices=SPAN_SIZES,
                    help="locus size at which --target-spanning applies")
    ap.add_argument("--out", default="downsample_plan.tsv")
    ap.add_argument("--drop", default="",
                    help="comma separated sample ids to exclude from the panel")
    a = ap.parse_args()

    with open(a.measured, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit("the measured table is empty")
    missing = [c for c in NEEDED if c not in rows[0]]
    if missing:
        sys.exit(f"the measured table lacks the columns {missing}")

    drop = {s.strip() for s in a.drop.split(",") if s.strip()}
    panel, dropped = [], []
    for r in rows:
        sid = r["label"]
        try:
            d = float(r["depth_x"])
        except ValueError:
            sys.exit(f"unreadable depth for {sid}: {r['depth_x']!r}")
        if d <= 0:
            sys.exit(f"non positive depth for {sid}: {d}")
        (dropped if sid in drop else panel).append((sid, d, r))

    if dropped:
        print(f"excluded on request: {' '.join(s for s, _, _ in dropped)}")
    if not panel:
        sys.exit("every sample was excluded")

    depths = [d for _, d, _ in panel]
    n50s = [float(r["read_n50"]) for _, _, r in panel]
    print(f"panel: {len(panel)} samples")
    print(f"  depth    median {statistics.median(depths):7.1f}x   "
          f"min {min(depths):6.1f}   max {max(depths):6.1f}")
    print(f"  read N50 median {statistics.median(n50s)/1000:7.1f} kb  "
          f"min {min(n50s)/1000:6.1f}   max {max(n50s)/1000:6.1f}")

    given = sum(x is not None and x is not False
                for x in (a.target, a.target_spanning, a.target_from_min or None))
    if given != 1:
        sys.exit("give exactly one of --target, --target-spanning or --target-from-min")

    span_col = f"span_{a.at_locus_bp // 1000}kb"
    basis = "depth_x"
    if a.target_spanning is not None:
        if span_col not in rows[0]:
            sys.exit(f"the measured table lacks the column {span_col}. Re-run the "
                     "measurement with the current scripts/measure_bam_lengths.sh, "
                     "which emits the spanning profile in the same pass.")
        basis = span_col
        target = a.target_spanning
        vals = []
        for sid, _, r in panel:
            try:
                v = float(r[span_col])
            except ValueError:
                sys.exit(f"unreadable {span_col} for {sid}: {r[span_col]!r}")
            if v <= 0:
                sys.exit(f"non positive {span_col} for {sid}: {v}. That sample cannot "
                         f"genotype a locus of {a.at_locus_bp} bp at all.")
            vals.append(v)
        print(f"  spanning reads at {a.at_locus_bp} bp, median "
              f"{statistics.median(vals):7.1f}   min {min(vals):6.1f}   max {max(vals):6.1f}")
        print(f"\ntarget: {target:.2f} reads spanning a {a.at_locus_bp} bp locus")
        print("  Equalising this rather than the raw depth is what makes the panel "
              "comparable,\n  but note it cannot equalise spanning depth at every locus "
              "size at once, since\n  the read length distributions differ.")
    elif a.target_from_min:
        target = min(depths)
        print(f"\ntarget taken from the shallowest sample: {target:.1f}x")
    else:
        target = a.target
        print(f"\ntarget given: {target:.1f}x")

    def value_of(r):
        return float(r[basis])

    below = [(s, value_of(r)) for s, _, r in panel if value_of(r) < target * 0.98]
    if below:
        print(f"\n{len(below)} samples sit below the target and cannot reach it by subsampling.")
        for s, d in sorted(below, key=lambda t: t[1])[:12]:
            print(f"    {s:<9} {d:8.2f}")
        print("  They are kept whole, so the panel is not uniform at the target. Either lower")
        print("  the target to the minimum with --target-from-min, or drop these samples.")

    out, fracs = [], []
    for sid, d, r in sorted(panel):
        f = min(1.0, target / value_of(r))
        note = "ok" if f < 0.98 else "kept_whole_at_or_below_target"
        out.append([sid, f"{d:.3f}", f"{f:.4f}", r["read_n50"],
                    r.get(span_col, ""), basis, f"{target:.4f}", note])
        fracs.append(f)

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "measured_depth_x", "subsample_fraction",
                    "read_n50", span_col, "target_basis", "target_value", "note"])
        w.writerows(out)

    kept = statistics.fmean(fracs)
    print(f"\nwritten {a.out}: {len(out)} samples")
    print(f"  fraction median {statistics.median(fracs):.3f}  "
          f"range {min(fracs):.3f} - {max(fracs):.3f}")
    print(f"  mean fraction retained {kept:.3f}, so the LAST stage falls to about "
          f"{kept*100:.0f}% of its native cost")
    print(f"  retained bases {sum(d*f for (_, d, _), f in zip(sorted(panel), fracs)):.0f}x "
          f"summed, against {sum(depths):.0f}x native")
    if basis != "depth_x":
        res = [float(r[span_col]) * f for (_, _, r), f in zip(sorted(panel), fracs)]
        print(f"  spanning reads at {a.at_locus_bp} bp after subsampling: "
              f"{min(res):.2f} to {max(res):.2f}, target {target:.2f}")


if __name__ == "__main__":
    main()
