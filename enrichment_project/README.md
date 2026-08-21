# Bulk RNA-seq — differential expression and enrichment

DESeq2 differential expression between breast cancer subtypes, followed by
over-representation analysis and GSEA. Python throughout (`pydeseq2`,
`gseapy`) rather than the usual R stack, so the whole thing runs in one
environment.

## Two scripts

| Script | Does | Writes |
| --- | --- | --- |
| `deseq2_analysis.py` | pre-filter → `DeseqDataSet(design="~Subtype")` → size factors → Wald test → LFC shrinkage | `deseq2_full_results.csv`, `overexpressed_genes.csv`, `downexpressed_genes.csv`, volcano and boxplot figures |
| `enrichment_analysis.py` | ORA on the up/down gene lists, GSEA on the full ranked result, PCA on normalised counts | `gsea_summary.png`, `enrichment_emt.png`, `pca_plot.png` |

The second reads the CSVs the first writes, so they run in order:

```bash
python deseq2_analysis.py
python enrichment_analysis.py     # or: sbatch enrichment_analysis.sh
```

## Inputs

Not versioned here — `dataset/` holds `STAR_counts.tsv`, `metadata.tsv` and
`gencode.v42.genes.tsv`.

Pre-filtering keeps genes with ≥ 10 counts in ≥ 3 samples. ORA queries
`GO_Molecular_Function_2023`, `GO_Cellular_Component_2023`,
`GO_Biological_Process_2023` and `KEGG_2021_Human` via Enrichr; GSEA is run on
genes ranked by the Wald statistic.

## Figures

| File | Shows |
| --- | --- |
| `volcano_(raw).png`, `volcano_(shrunk).png` | DE before and after LFC shrinkage |
| `raw_vs_shrunk_lfc.png` | what shrinkage does to low-count genes |
| `top4_degs.png`, `esr1_boxplot.png` | top hits, and ESR1 as the subtype sanity check |
| `pca_plot.png` | samples on normalised counts |
| `gsea_summary.png`, `enrichment_emt.png` | GSEA summary and the EMT gene set |

Both volcano versions are kept deliberately: the raw plot is what makes
shrinkage look like it is hiding signal, and the pair together shows it is
removing noisy low-count log-fold-changes instead.

Note: comments in both scripts are in Russian.
