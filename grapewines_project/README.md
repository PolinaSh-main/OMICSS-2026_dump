# Grapevine population genomics — OMICSS 2026, GP3

Population structure and differentiation of Caucasian grapevine
(*Vitis vinifera*) from SNP data: 412 samples, ~815 k biallelic SNPs after
filtering, 19 chromosomes.

Course instructions live in the
[OMICSS repo](https://github.com/abi-am/omicss26/tree/main/Group%20projects/GP3).
This checkout turns them into pipelines that can be re-run.

## Layout

```
grapewines_project/
├── conf/base.config          shared: SLURM executor, reports, common params
├── reference/                small inputs that belong in git
│   ├── cauc_filtered.final.fam
│   └── cauc_grape_metadata.csv
├── pipelines/
│   ├── 00_preprocess/        VCF filtering, VCF -> PLINK
│   ├── 01_admixture/         CV curve, sample orderings, barplot grid
│   ├── 02_pca/               PLINK PCA + scatter plots per metadata column
│   ├── 03_fst/               pairwise FST between ADMIXTURE groups
│   ├── 04_manhattan_annotation/  windowed FST -> outlier windows -> genes
│   └── 05_tree/              neighbour-joining tree from 1-IBS distances
└── results/                  published outputs (small ones are versioned)
    ├── admixture/{cv,orders,orderings,barplots}
    ├── pca/{global,by_metadata,association}
    ├── fst/K<K>/
    ├── manhattan/<comparison>/
    └── tree/
```

Every pipeline is launched from its own directory, so `projectDir` resolves
correctly and `bin/` lands on `PATH`:

```bash
cd pipelines/03_fst
nextflow run . -profile slurm --k 7
```

Nothing is hardcoded to one user's home any more: `conf/base.config` derives
everything from `params.project_root`, which defaults to two levels above the
pipeline directory. Override it with `--project_root` if needed.

## Data

| What | Where | In git? |
|---|---|---|
| filtered VCF + index | `/mnt/nas1/proj/omicss26/gp3/data/vcf/cauc_filtered.final.vcf.gz{,.tbi}` | no (shared, read-only) |
| PLINK triple | `<...>/data/plink/cauc_filtered.final.{bed,bim,fam}` | `.fam` only |
| metadata | `reference/cauc_grape_metadata.csv` | yes |
| ADMIXTURE `.Q` / `.P` / logs | user scratch, see `pipelines/01_admixture/nextflow.config` | no |
| gene annotation reference | `/mnt/nas1/proj/omicss26/gp3/gene_annotation/reference/PN40024.v4.1.REF.{gff3,b2g.tsv}` | no |

## Analysis order

```
00_preprocess   VCF -> filtered VCF -> PLINK bed/bim/fam
      |
01_admixture    ADMIXTURE K=2..10 -> CV curve picks K -> barplots
      |
02_pca          independent cross-check of the same structure
      |
03_fst          assign samples to K groups -> pairwise FST -> heatmap
      |
04_manhattan_annotation
                per-SNP FST -> 50 kb windows -> top 1% -> genes in those windows

05_tree         independent of all of the above: NJ tree on 1-IBS distances
```

## Tooling

Everything is Python — pandas, numpy, matplotlib — plus the usual
command-line tools (vcftools, PLINK, ADMIXTURE) driven by Nextflow. The
course materials write the FST heatmap and the Manhattan plot in R;
those are Python here instead, so the project needs one environment
rather than two and the plotting code can share the palette below with
the ADMIXTURE figures. The outputs are the same files.

## Colour palette

Shared between ADMIXTURE barplots, metadata barplots and FST figures so that
component *n* is the same colour everywhere:

| | | | | | | | |
|---|---|---|---|---|---|---|---|
| K1 `#00699A` | K2 `#FD702B` | K3 `#06402B` | K4 `#FDC700` | K5 `#AE0039` | K6 `#6A005C` | K7 `#08BDBD` | K8 `#F21B3F` |

## Notes

- `pipelines/00_preprocess/bin/filter_snp_*.sh` and `vcf_to_plink.sh` are the
  instructors' original scripts, kept for provenance. They still contain
  absolute paths from the shared project folder.
- `pipelines/02_pca/bin/run_pca_plots.sh` calls `pca_plot.py` without
  `--color-column`, which that script requires — it fails as written. Use
  `run_all_pca_plots.sh`, which passes the column.
- `pipelines/02_pca/bin/pc_metadata_heatmap.py` asks which metadata
  variable explains which PC — one-way ANOVA per PC × variable, cells
  coloured by adjusted R², stars from the BH-adjusted p-value. It is a
  port of a `heatmap.r` written in RStudio; the R on this cluster has no
  ggplot2 and no png device, so it could not have run here. The port was
  checked against R's own `lm`/`anova` on this data and agrees to 1e-10.
  Defaults are sized for a slide — PC1–PC5 across eight variables, large
  type; pass `--font-scale 0.6` for a figure that goes in a document.

  ```bash
  cd pipelines/02_pca
  python3 bin/pc_metadata_heatmap.py \
      --eigenvec ../../results/pca/global/cauc_pca.eigenvec \
      --metadata ../../reference/cauc_grape_metadata.csv \
      --outdir   ../../results/pca/association
  ```
- `01_admixture` was moved from `admixture/` into `pipelines/01_admixture/`
  and its config was made relocatable; the results in `results/admixture/`
  come from the run performed before that move.
