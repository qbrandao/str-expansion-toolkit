# str-expansion-toolkit

CLI pour la détection et l'analyse d'expansions de répétitions courtes (STR)
à partir de données long-read ONT, en combinant plusieurs outils orthogonaux :
**VAMOS**, **tandem-genotypes**, **LongTR** (par défaut), et **TRGT** (en
option, voir note ci-dessous).

## ⚠️ Note importante sur TRGT et les données ONT

TRGT est conçu pour des reads PacBio HiFi et **n'a pas de support officiel
pour les données ONT** (cf. Aliyev et al. 2026, bioRxiv, qui exclut
explicitement TRGT d'un benchmark d'outils de génotypage STR sur données
Nanopore pour cette raison — il n'accepte pas les données ONT en usage
prévu). En conséquence, **TRGT n'est pas lancé par défaut** par `detect`. Il
reste disponible via `--tools ... trgt ...` pour comparaison/exploration,
mais tout résultat produit ainsi doit être documenté comme un usage hors
cadre officiel dans toute publication qui en découlerait.

Le 3e outil par défaut est **LongTR**, une adaptation de HipSTR conçue pour
les reads longs PacBio HiFi **et** ONT.

## Installation

```bash
git clone https://github.com/qbrandao/str-expansion-toolkit.git
cd str-expansion-toolkit
pip install -e .
```

## Configuration

Tous les chemins spécifiques à l'installation (génome de référence, catalogue
de motifs VAMOS, index minimap2 pour TRGT/LongTR, index LAST pour
tandem-genotypes, modèle clair3, noms des environnements micromamba) sont
centralisés dans `config.yaml`, à partir du modèle `config.example.yaml` :

```bash
cp config.example.yaml config.yaml
# éditer config.yaml avec tes propres chemins
```

## Utilisation

### 1) Détecter les expansions chez un ou plusieurs patients

```bash
# Un seul patient (outils par défaut : vamos, tandem-genotypes, longtr)
str-toolkit detect \
  --sample patient01 \
  --bam patient01.sorted.bam \
  --fastq patient01_Guppy_4.0.11_prom.merged.fastq.gz \
  --config config.yaml \
  -o results/patients/

# Une liste de patients (TSV: sample_id, bam_path, fastq_path)
str-toolkit detect \
  --samples-list patients.tsv \
  --config config.yaml \
  -o results/patients/

# Inclure TRGT explicitement (hors cadre officiel sur ONT, à ne faire
# qu'en connaissance de cause -- voir note en haut de ce fichier)
str-toolkit detect \
  --sample patient01 --bam ... --fastq ... --config config.yaml \
  --tools vamos tandem-genotypes longtr trgt \
  -o results/patients/
```

- `--bam` : BAM déjà aligné, utilisé par VAMOS (clair3 + whatshap + vamos --contig).
- `--fastq` : fastq(.gz) brut fusionné, utilisé par TRGT/LongTR (alignement
  minimap2, préréglage `map-ont`) et tandem-genotypes (last-train/lastal).
  TRGT et LongTR réutilisent automatiquement le même `.sorted.bam` s'ils
  tournent tous les deux dans le même run (pas de double alignement).

Lance les outils sélectionnés pour chaque échantillon (avec reprise sur
erreur : les étapes dont la sortie existe déjà sont sautées), puis fusionne
les sorties en un VCF récapitulatif : `results/patients/<sample_id>/<sample_id>.merged.vcf`.

**Fusion (`str_toolkit/merge.py`) :** les outils n'ancrent pas leurs
coordonnées de la même façon pour un même locus biologique. La fusion
regroupe donc les appels par **intervalle tolérant** (`window`, 25 bp par
défaut) **et motif canonique** (rotation circulaire : `AAAG` ≡ `AAGA` ≡
`GAAA`), plutôt que par égalité stricte de position. Les tailles restent
**séparées par source** dans le VCF (`SIZES=vamos_hap1:42,longtr_allele1:12,...`)
car leurs unités ne sont pas comparables :
- VAMOS : longueur en unités de motif
- TRGT : longueur en bp (allèle absolu)
- LongTR : différence en bp par rapport à la référence (delta, pas une longueur absolue)
- tandem-genotypes : longueur en bp dérivée du clustering de longueurs read-level

### 2) Construire le référentiel de contrôles

```bash
str-toolkit build-controls \
  --controls-dir results/controls/ \
  -o results/controls.json
```

Lit le VCF fusionné de chaque contrôle
(`{controls-dir}/{sample_id}/{sample_id}.merged.vcf`, produit par `detect`) et
construit, pour chaque locus, la taille maximale observée **séparément par
outil** :

```json
{
  "chr1_12345_AAAG": {
    "chrom": "chr1", "pos": 12345, "motif": "AAAG",
    "tools": {
      "vamos": {"max_size": 42, "n_observed": 87},
      "tandem-genotypes": {"max_size": 3, "n_observed": 85},
      "longtr": {"max_size": 10, "n_observed": 88}
    }
  }
}
```

Optionnel : `--samples-list` pour restreindre aux échantillons listés (sinon
tous les sous-dossiers de `--controls-dir` sont utilisés).

### 3) Comparer les patients aux contrôles

```bash
str-toolkit compare \
  --patients-dir results/patients/ \
  --controls-json results/controls.json \
  --genes-bed genes.bed.gz \
  --exons-bed MANE_Select_exons.bed.gz \
  -x 0 \
  -o report.tsv
```

Pour chaque locus patient présent dans le référentiel de contrôles, calcule
un diff **par outil disponible** (`vamos_diff`, `trgt_diff`,
`tandem_genotypes_diff`, `longtr_diff`), garde la ligne si au moins un outil
dépasse le seuil (`-x/--threshold`), et trie décroissant sur `max_diff` (le
plus grand diff parmi les outils disponibles pour ce locus).

Colonnes de sortie :
`patient_id, chrom, pos, motif, gene, feature, vamos_size, vamos_control_max,
vamos_diff, trgt_size, trgt_control_max, trgt_diff, tandem_genotypes_size,
tandem_genotypes_control_max, tandem_genotypes_diff, longtr_size,
longtr_control_max, longtr_diff, n_tools_expanded, max_diff`

`n_tools_expanded` compte combien d'outils orthogonaux confirment l'expansion
à ce locus — un signal de confiance utile (une expansion vue par plusieurs
outils indépendants est plus fiable qu'une vue par un seul).

Options : `-t/--triplet-only` (motifs ≥ 3 pb uniquement), `--patients`
(restreindre à une liste d'identifiants), `--format csv`.

## Format des fichiers d'échantillons

Fichier TSV avec en-tête, utilisé par `--samples-list` de `detect` (`bam_path`
et/ou `fastq_path` selon les outils utilisés) :

```
sample_id	bam_path	fastq_path
patient01	/data/bam/patient01.sorted.bam	/data/fastq/patient01.merged.fastq.gz
patient02	/data/bam/patient02.sorted.bam	/data/fastq/patient02.merged.fastq.gz
```

## Statut du projet

Fonctionnel de bout en bout : `detect` (VAMOS/tandem-genotypes/LongTR par
défaut, TRGT en opt-in + fusion multi-outils), `build-controls` et `compare`
(registre et diffs par outil, annotation gène/feature).

Formats de sortie **confirmés** sur des fichiers réels :
- TRGT : `INFO/MOTIFS`, `FORMAT/AL` (longueurs d'allèles en bp).
- tandem-genotypes : TSV à 8 colonnes (`chrom, start, end, motif, (ignorée), .,
  longueurs_par_read, .`). La colonne des longueurs par read est séparée en
  2 groupes (allèle court/long) au niveau du plus grand écart, dont on
  prend la médiane comme taille d'allèle (bp). Le bed TRF listant souvent
  plusieurs motifs candidats qui se chevauchent pour un même locus, ces
  candidats sont dédoublonnés en gardant celui couvert par le plus de reads.
- LongTR : `INFO/MOTIF`, `FORMAT/GB` (différence en bp vs référence, par
  allèle) — d'après le README officiel gymrek-lab/LongTR.

Les environnements micromamba référencés dans `config.yaml` (clair3,
whatshap-env, vamos, trgt, longtr, last_env, tandem-env) doivent déjà exister
sur la machine d'exécution.

## Licence

MIT
