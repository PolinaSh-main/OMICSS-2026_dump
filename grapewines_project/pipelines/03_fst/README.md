# 03_fst — differentiation between ADMIXTURE groups

ADMIXTURE says how much ancestry each sample draws from each cluster.
It does not say how far apart those clusters are. FST does.

```
.Q + .fam  ->  sample lists   ->  pairwise FST  ->  summary  ->  heatmap
                    |
                metadata  ->  proposed group names
```

## Run

```bash
cd pipelines/03_fst

nextflow run . --k 7                 # all 21 comparisons
nextflow run . --k 3                 # the instructors' worked example
nextflow run . --k 7 --pairs K1:K5   # one comparison
```

Everything runs through SLURM. The 21 comparisons at K=7 are independent
and start as 21 separate jobs, so the wall time is one vcftools pass over
the VCF rather than twenty-one of them in sequence.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--k` | 7 | which ADMIXTURE run to use |
| `--admixture_dir` | user scratch | directory holding `cauc_filtered.final.<K>.Q` |
| `--min_q` | 0.75 | ancestry proportion required to assign a sample |
| `--min_group_size` | 4 | groups smaller than this are not compared |
| `--pairs` | `all` | `all`, or `K1:K5`, or `K1:K5,K2:K3` |
| `--metric` | `mean_fst` | column the heatmap shows |

## Output

Written to `results/fst/K<K>/`, in the layout the course instructions ask
for:

```
sample_lists/
    K1_samples.txt ... K7_samples.txt   sample IDs per group
    admixed_samples.txt                 below the 0.75 threshold
    sample_q_values_K7.tsv              every sample + its Q values
    group_sizes_K7.tsv                  deliverable: group sizes
    group_interpretation_K7.tsv         deliverable: proposed names
    group_composition_K7.tsv            country/utilization breakdown
    group_labels_K7.tsv                 axis labels for the heatmap
pairwise_fst/
    K1_vs_K2.weir.fst  ...              per-SNP FST (git-ignored, large)
    K1_vs_K2.log       ...              vcftools logs
    pairwise_fst_summary_K7.tsv         deliverable: one row per pair
plots/
    pairwise_fst_heatmap_K7.pdf / .png  deliverable
    pairwise_fst_matrix_K7.csv
```

## Two FST numbers, not one

`pairwise_fst_summary_K<K>.tsv` reports both averages vcftools makes
available:

- `mean_fst` — plain average over SNPs. This is what the course
  instructions ask for, and what the heatmap shows by default.
- `weighted_fst` — the ratio-of-averages estimator. Less swayed by the
  very many near-zero SNPs, and the one usually quoted in papers.

They differ by roughly a third on this data. Say which one a figure
shows.

## Reading the heatmap

The diagonal is zero by construction and is left blank. The colour ramp
is stretched across the range of the values that are actually present,
not from zero — otherwise, when all the pairwise values sit close
together, every cell comes out the same shade. **The colours are
relative to this figure.** Every cell carries its number, and the colour
bar is annotated, so read those against the usual guide (≈0.05 low,
≈0.15 moderate, >0.25 strong) rather than judging by shade alone.

## Naming the groups

`describe_groups.py` proposes a name for each component from the
metadata — dominant country, plus wild or cultivated inferred from
`Genetic background` (ssp. *sylvestris* = wild, ssp. *vinifera* =
cultivated) with `Utilization` as a fallback. Two components that land on
the same label get numbered, as in the instructors' "Armenian wild group
1 / 2".

The `evidence` column carries the numbers behind each name, and
`homogeneous` flags groups where no single country reaches 60%. **Read
these before using the names in a report** — the script is a first pass
over the metadata, not a conclusion.
