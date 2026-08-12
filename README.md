# OMICSS-2026

Course work from the OMICSS 2026 summer school. Four independent projects,
each self-contained.

| Directory | Topic |
|---|---|
| [`grapewines_project/`](grapewines_project/) | **Group project 3** — population genomics of Caucasian grapevine (PLINK / ADMIXTURE / PCA / FST / gene annotation). This is the active one. |
| [`tmf_project/`](tmf_project/) | WES variant calling practice (BWA → GATK) |
| [`single_cell_pipeline/`](single_cell_pipeline/) | scRNA-seq practice (QC → doublets → Harmony → UMAP) |
| [`enrichment_project/`](enrichment_project/) | Bulk RNA-seq: DESeq2 + GSEA practice |

## What is in git and what is not

Git holds **code, configs and small text results** — nothing that a pipeline
can regenerate.

Kept: `*.nf`, `*.config`, scripts, `.fam`, `.Q`, summary tables, final figures.

Ignored (see [`.gitignore`](.gitignore)): `work/`, `.nextflow/`, cluster logs,
Nextflow HTML reports, virtualenvs, and all bulk genomics data — FASTQ, BAM,
VCF, PLINK `.bed`/`.bim`, ADMIXTURE `.P`, per-SNP `.weir.fst`, `.h5ad`.

These files still live on the server; they are simply not versioned.

## Working on the server

The repository is checked out at

```
/mnt/nas0/user/polina.shevyakova/
```

Course-wide read-only data provided by the instructors is at

```
/mnt/nas1/proj/omicss26/gp3/
```
