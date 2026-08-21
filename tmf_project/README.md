# WES variant calling — BWA → GATK → ANNOVAR

Whole-exome sequencing of two samples (`wes46`, `wes78`), taken from aligned
BAMs through to an annotated, human-readable variant table for one gene of
interest (`MEFV`, familial Mediterranean fever).

## Pipeline

```
BAM ──▶ SORT_BAM ──▶ MARK_DUPLICATES ──▶ HAPLOTYPECALLER (per sample, gVCF)
                                              │
                                    COMBINE_GVCFS ──▶ GENOTYPE_GVCFS
                                              │
                                      FILTER_VARIANTS  (SNPs / indels split)
                                              │
                                         ANNOTATE      (ANNOVAR, hg38)
                                         │        │
                                  ADD_RSIDS    BEAUTIFY_TABLE
```

| Stage | Tool | Output |
| --- | --- | --- |
| `SORT_BAM`, `MARK_DUPLICATES` | GATK | coordinate-sorted, deduplicated BAM |
| `HAPLOTYPECALLER` | GATK | per-sample `.g.vcf.gz` |
| `COMBINE_GVCFS`, `GENOTYPE_GVCFS` | GATK | joint-genotyped multi-sample VCF |
| `FILTER_VARIANTS` | GATK `VariantFiltration` | `filtered_snps.vcf`, `filtered_indels.vcf` |
| `ANNOTATE` | ANNOVAR (`refGeneWithVer`, `avsnp150`) | `*_multianno.{txt,vcf}` |
| `ADD_RSIDS` | bcftools | VCFs with `ID` set from `avsnp150` |
| `BEAUTIFY_TABLE` | shell | `readable_{snps,indels}.txt`, `MEFV_report.txt` |

Every process wraps its command in a per-sample log under
`results/<runName>/logs/`, so a failure is traceable to one stage and one
sample without re-running anything.

`HAPLOTYPECALLER` takes a `fallback` gVCF path per sample: exome-scale calling
is the slowest stage by far, so a pre-computed gVCF can be substituted when one
is available rather than recomputing it.

## Inputs

| What | Where | In git? |
| --- | --- | --- |
| aligned BAMs | `${params.project}/bam/<sample>.bam` | no |
| hg38 reference | `params.ref` | no |
| ANNOVAR + `humandb` | `params.annovar` | no |
| GATK 4.2.6.1 | `params.gatk` | no |

Paths in `nextflow.config` are placeholders (`/path/to/...`) — point them at
your own installation, or override on the command line.

## Running

```bash
cd variant_calling
nextflow run . --project /path/to/workdir --ref /path/to/hg38.fa
```

The config targets a SLURM executor with per-process cpu/memory/time requests.
Drop `executor = 'slurm'` from `nextflow.config` to run locally.

## Also here

`scr/` holds the standalone shell steps used before the pipeline was written —
`bwa_index.sh`, `alignment.sh` (BWA-MEM), `run_fastqc.sh`, `run_multiqc.sh`,
`gatk_sort.sh`. FastQC/MultiQC reports for the four read files are kept under
`data/fastq/qc_result/`.
