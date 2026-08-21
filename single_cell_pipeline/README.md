# scRNA-seq — QC → doublet removal → integration

Two 10x conditions (`control`, `stim`) taken from raw count matrices to a
Harmony-integrated embedding, as a Nextflow DSL2 workflow with one Python
script per stage.

## Pipeline

```
control/ ──┐
           ├──▶ QC ──▶ DOUBLETS ──┐
stim/    ──┘   (per condition)    ├──▶ MERGE ──▶ NORMALIZE ──▶ HARMONY
                                  ┘
```

| Stage | Script | Does |
| --- | --- | --- |
| `QC` | `scripts/qc.py` | filter cells on gene count and mitochondrial fraction |
| `DOUBLETS` | `scripts/doublet_removal.py` | per-condition doublet detection and removal |
| `MERGE` | `scripts/merge.py` | concatenate conditions into one AnnData |
| `NORMALIZE` | `scripts/normalize.py` | `normalize_total` → `log1p` → HVG → scale → PCA |
| `HARMONY` | `scripts/harmony.py` | batch-correct the PCA embedding across conditions |

`scripts/run_umap.py` and `scripts/plot_umap.py` produce the UMAP from the
integrated embedding; `scripts/visualize.py` holds the shared plotting code.

QC and integration are separated on purpose: doublet detection runs *per
condition*, before merging, because doublet rates are a property of a single
10x run and estimating them on a merged object conflates two libraries.

## Parameters

Defaults in `main.nf`, all overridable on the command line:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--min_genes` | 200 | drop cells below this many detected genes |
| `--max_genes` | 2500 | drop suspected multiplets above this |
| `--max_mt` | 5.0 | max % mitochondrial counts |
| `--expected_doublet_rate` | 0.06 | prior doublet rate |
| `--n_top_genes` | 2000 | highly variable genes kept |
| `--n_pcs` | 30 | principal components into Harmony |

## Inputs

10x matrices, one directory per condition:

```
data/control/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
data/stim/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
```

`barcodes.tsv.gz` and `features.tsv.gz` are small and versioned here; the
`matrix.mtx.gz` files are not — they are ~27 MB each and regenerable from the
source dataset.

## Running

```bash
python -m venv .venv && .venv/bin/pip install scanpy anndata harmonypy
cd single_cell_pipeline
nextflow run .
```

The processes call `${projectDir}/.venv/bin/python` directly rather than
relying on an activated environment, so the interpreter is the same whichever
node SLURM picks. Create that venv in the pipeline directory before running.
Drop `executor = 'slurm'` from `nextflow.config` to run locally.
