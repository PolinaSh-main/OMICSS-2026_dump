# 04_manhattan_annotation — from FST peaks to candidate genes

```
*.weir.fst  ->  50 kb windows  ->  Manhattan plot + top 1% windows
            ->  merged regions ->  genes in/near them  ->  report table
```

Input is whatever `03_fst` already produced. Nothing is recomputed from
the VCF.

## Run

```bash
cd pipelines/04_manhattan_annotation

nextflow run . --k 7                          # every comparison
nextflow run . --k 7 --comparisons K1_vs_K3   # just one
```

## Why windows and not SNPs

A single SNP with high FST means very little — with ~815 k SNPs, the
extreme tail is full of noise. Averaging over a 50 kb window and
requiring at least 20 SNPs in it turns the tail into something that
reflects a genuine stretch of differentiated genome. **One point on the
Manhattan plot is one window.**

The top 1% of windows by mean FST are marked. Adjacent outlier windows
are then merged: one sweep usually shows up as a run of high windows, and
reporting each of them separately would count the same locus several
times.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--k` | 7 | which `results/fst/K<K>/` to read |
| `--comparisons` | `all` | `all`, or `K1_vs_K3`, or a comma-separated list |
| `--window` | 50000 | window size, bp |
| `--min_snps` | 20 | windows with fewer SNPs are dropped |
| `--top_fraction` | 0.01 | fraction of windows called outlying |
| `--max_gap` | 50000 | outlier windows this close merge into one region |
| `--flank` | 10000 | also report genes this far outside a region |
| `--top_regions` | 5 | regions in the final table |
| `--genes_per_region` | 2 | genes per region in the final table |

## Output

Per comparison, in `results/manhattan/K<K>/<comparison>/`:

| File | What |
|---|---|
| `*_windows.tsv` | every window that passed `--min_snps` |
| `*_outlier_windows.tsv` | the top 1% |
| `*_manhattan.pdf` / `.png` | the plot |
| `*_candidate_regions.tsv` | merged regions, with `mean_window_fst` and `max_window_fst` |
| `*_candidate_genes.tsv` | every gene in or near a region, with `distance_to_region_bp` |
| `*_selected_regions.tsv` | **the report table** |
| `*_low_information_genes.tsv` | hypothetical / unnamed / transposon genes, set aside |

## How genes are picked

Taking the first rows of `*_candidate_genes.tsv` would give the same
locus several times over, since a region holds many genes. Instead:

1. distinct regions are ranked by `max_window_fst`, then `mean_window_fst`;
2. the top five are kept;
3. inside each, genes that overlap the region (`distance_to_region_bp = 0`)
   come first;
4. among those, an annotation matching disease resistance, stress
   response, hormone signalling, berry metabolism or development wins,
   and one or two are kept;
5. hypothetical, uncharacterised and transposon-related genes go to
   `*_low_information_genes.tsv`, and are only promoted into the report
   when a region has nothing better.

The `theme` and `note` columns record why each gene was chosen. The
themes are keyword matches on the Blast2GO description — a shortlist to
read, not a functional claim.

## Chromosome naming

The FST files come from the VCF and use `1`..`19`. The PN40024 reference
GFF3 uses `chr01`..`chr19`. Both sides are reduced to their digits before
matching, and the pipeline stops with an explicit error if that still
leaves no shared chromosome — rather than quietly reporting zero genes.
