# Scripts SLURM de test

À lancer depuis la racine du repo (`sbatch scripts/00_unit_tests.sbatch`, etc.),
ou tout enchaîner d'un coup avec `./scripts/submit_full_pipeline.sh`.

## Prérequis (à placer à la racine du repo avant de lancer)

- `config.yaml` (copié depuis `config.example.yaml`, chemins adaptés)
- `controls.tsv`, `patients.tsv` : TSV avec en-tête `sample_id\tbam_path\tfastq_path`
- `genes.bed.gz`, `MANE_Select_exons.bed.gz` : utilisés par `compare`
- `longtr.regions_bed` (config.yaml) : voir `./scripts/prepare_longtr_bed.sh` ci-dessous
  pour le générer à partir du catalogue Adotto (nécessite `gawk`)
- Adapter `<ta_partition>` dans chaque `.sbatch` (partition SLURM de ton cluster)
- `mkdir -p logs` (fait automatiquement par les scripts, mais le dossier
  doit être accessible en écriture depuis les nœuds de calcul)

## Préparer le catalogue de régions LongTR

```bash
./scripts/prepare_longtr_bed.sh /mnt/references
```

Télécharge le catalogue Adotto (hg38) et le convertit au format attendu par
LongTR (`chrom, start 1-based, end, motif[,motif2], nom`). Nécessite `gawk`
(le script vérifie sa présence et s'arrête sinon).

## Ordre recommandé (première fois, pas à pas)

1. `00_unit_tests.sbatch` — valide l'installation du package (aucun outil bioinfo lancé)
2. `01_smoke_test_single_tool.sbatch` — valide `config.yaml`/les chemins avec TRGT seul (rapide)
3. `02_detect_single_patient.sbatch` — valide le pipeline complet (3 outils) sur 1 patient
4. `03_detect_controls_array.sbatch` — détection sur toute la cohorte contrôle (job array,
   adapter `--array=1-N%4` où N = nombre de lignes de `controls.tsv`)
5. `04_detect_patients_array.sbatch` — détection sur tous les patients (job array, idem)
6. `05_build_controls.sbatch` — après la fin complète de l'étape 4 (registre JSON)
7. `06_compare.sbatch` — après la fin des étapes 5 et 6 (rapport final)

## Tout enchaîner automatiquement

```bash
./scripts/submit_full_pipeline.sh
```

Soumet 00 → 01 → (03 array + 04 array en parallèle) → 05 → 06, avec les
dépendances `--dependency=afterok:...` calculées automatiquement (nombre de
lignes de `controls.tsv`/`patients.tsv` compté dynamiquement). Suivi avec
`squeue -u $USER`.
