# 05_tree — neighbour-joining tree

```
PLINK bed/bim/fam
     -> LD pruning (r2 < 0.2)
     -> 1-IBS distance matrix + genotype export
     -> NJ tree, bootstrapped, midpoint rooted
     -> circular figures coloured by metadata
```

A third, independent view of the same structure that ADMIXTURE and PCA
describe. ADMIXTURE assumes K ancestral populations; PCA assumes the
structure is linear; NJ assumes neither, so agreement between the three
is worth more than any of them alone.

## Run

```bash
cd pipelines/05_tree

nextflow run .

# colour by ADMIXTURE group instead of metadata
nextflow run . \
    --color_by assignment \
    --groups ../../results/fst/K7/sample_lists/sample_q_values_K7.tsv

# quick look, no bootstrap
nextflow run . --bootstrap 0
```

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--bfile` | shared course PLINK set | prefix, without `.bed/.bim/.fam` |
| `--chr_set` | 19 | grapevine chromosome count |
| `--ld_r2` | 0.2 | LD pruning threshold |
| `--bootstrap` | 100 | replicates; 0 to skip |
| `--root` | `midpoint` | `midpoint`, `outgroup`, `none` |
| `--outgroup` | — | sample ID, when `--root outgroup` |
| `--color_by` | 3 metadata columns | comma-separated; `assignment` uses the ADMIXTURE group |
| `--groups` | — | `sample_q_values_K<K>.tsv`, for `assignment` |
| `--layout` | `circular` | `circular` or `rectangular` |
| `--scale` | `linear,cladogram` | comma-separated; one figure each |
| `--show_labels` | false | tip names; unreadable above ~80 samples |

## Output

In `results/tree/`:

| File | What |
|---|---|
| `tree.nwk` | Newick, bootstrap support on internal nodes |
| `distance_matrix.tsv` | the 1-IBS matrix |
| `tree_tips.tsv` | tip order |
| `plots/tree_<column>.pdf` / `.png` | branch lengths to scale |
| `plots/tree_<column>_cladogram.pdf` / `.png` | topology only |

## Why LD pruning first

Linked SNPs carry the same signal several times. Left in, they stretch
branch lengths in whichever regions happen to be dense, and — worse —
make bootstrap support look far better than it is, because a resample
that picks one SNP of a linked block effectively picks all of them. The
pruning step is what makes the support values mean anything.

## Why the linear tree looks like a starburst

Each sample has its own long private branch, because an individual
carries a lot of variation nobody else does. The splits *between* groups
are short by comparison. So on a true-to-scale drawing every tip sits at
roughly the same radius and the interesting structure is squeezed into a
small disc in the middle.

That figure is honest and should go in the appendix. The cladogram is
the one to read the topology from — its title says branch lengths are
not to scale, and it must stay saying that.

## Rooting

There is no outgroup in this dataset; every sample is *Vitis vinifera*.
The tree is therefore **midpoint rooted**, which places the root halfway
along the longest tip-to-tip path. That is a drawing convention chosen
so the fan has a centre — it is not a claim about which lineage is
ancestral. Do not describe any clade as "basal" on the strength of it.

## Distances, and why they are computed twice

PLINK produces the `1-IBS` matrix, and `build_tree.py` computes the same
matrix itself from the genotypes. The two are compared, and the pipeline
stops if they differ by more than 1e-4.

This is not redundancy for its own sake: the bootstrap has to rebuild
distances from resampled loci, which PLINK will not do, so that
estimator has to exist in Python anyway. Checking it against PLINK on
the full data is what makes it trustworthy.

## Reading the support values

Most internal nodes will have low support, and that is the correct
answer, not a bug. With 412 individuals, the majority of nodes describe
which of two nearly identical vines happens to pair with which — there
is no signal there to support. What matters is the support on the deep
splits that separate the groups. Quote those, and say plainly that the
within-group topology is unresolved.
