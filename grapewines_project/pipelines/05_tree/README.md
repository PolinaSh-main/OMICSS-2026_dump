# 05_tree — rooted maximum-likelihood phylogeny

Genome-wide SNP tree of the 412 Caucasian accessions with **ZZ01**
(*Vitis rotundifolia*) as an external outgroup, built with SNPhylo and
annotated with the K = 7 ADMIXTURE result.

```
cd pipelines/05_tree
nextflow run .
```

The maximum-likelihood stage takes hours. To start from a tree that has
already been built:

```
nextflow run . \
    --ml_tree /path/to/cauc_rooted.ml.tree \
    --phylip  /path/to/cauc_rooted.phylip.txt
```

## Steps

| process | what it does |
|---|---|
| `PREFLIGHT` | 413 samples, ZZ01 exactly once, 119 195 variants, no duplicate names |
| `SNPHYLO_ML` | LD pruning and SNP sequence via SNPRelate, then PHYLIP `dnaml` rooted on ZZ01 |
| `BOOTSTRAP` | 100 replicates with `phangorn::bootstrap.pml`, NNI optimisation |
| `VALIDATE` | both trees non-empty, 413 unique tips, ZZ01 present, support labels exist |
| `ANNOTATE` | tip table (K group / Admixed / Outgroup, max Q) and group descriptions |
| `PLOT` | fan and rectangular figures, PDF and PNG |

## Fixed settings

MAF ≥ 0.10, missing rate ≤ 0.10, LD 0.10, chromosomes 1–19, outgroup
ZZ01, 100 bootstrap replicates, low-depth screen skipped (`-r`). These
come from the task and are not to be changed without the mentor
agreeing; an approved alternative is a separate, labelled sensitivity
run — see `05b_nj_tree` for one such.

## Four things that are easy to get wrong here

**The SNPhylo path in the task brief is not readable.**
`/mnt/nas0/proj/vine/user_projects/armaria/snphylo_prog/` sits under a
directory owned by group `vine`. This pipeline uses the copy inside the
course project, `phylogeny/software/snphylo/`, which is the one the
group's earlier runs used and whose `snphylo.vcf.sh` has `BASE_DIR`
hard-coded to that location.

**LD 0.10 is not a stylistic choice.** SNPhylo refuses to build a tree
from an alignment longer than 50 000 sites. The group's earlier attempt
used 0.18, kept 55 275 sites and died on exactly that check. At 0.10
this input keeps **12 848**.

**SNPhylo's own bootstrap script oversubscribes the node.**
`determine_bs_tree.R` calls `bootstrap.pml(..., multicore = TRUE)`
without `mc.cores`, so phangorn falls back to `detectCores()` — 64 on
these machines, whatever SLURM allocated. A job asking for 5 CPUs forks
64 workers, each with its own copy of the alignment. That is how the
earlier run was killed. `bin/bootstrap_tree.R` passes the count
explicitly.

**K numbers are not group names.** ADMIXTURE orders its columns
arbitrarily, so `K3` means whatever this particular run made it mean.
The legend text is derived from the metadata in `annotate_tips.py`, not
copied from any table of K numbers belonging to a different run. On this
run the groups come out as:

| group | n | what it is |
|---|---|---|
| K1 | 15 | Armenian wild, group 1 |
| K2 | 22 | Georgian cultivated |
| K3 | 13 | Armenian–Azerbaijani cultivated, group 1 |
| K4 | 12 | Armenian–Azerbaijani cultivated, group 2 |
| K5 | 32 | Turkish wild |
| K6 | 49 | Turkish cultivated |
| K7 | 17 | Armenian wild, group 2 |
| Admixed | 252 | max Q below 0.75 |

**252 of 412 accessions — 61% — fall below the threshold.** They stay on
the tree in grey, and where they sit is worth reading, but only 160
samples carry a group colour.

## Figures

`tree_fan.*` shows all 413 accessions at once; tip labels are off,
because 413 names around a circle hide the pattern the figure exists to
show. `tree_rectangular.*` labels every tip, keeps branch lengths, and
prints bootstrap values of 70 or more.

The fan compresses the root stem. Rooting halves the long ZZ01 branch
and puts one half under each side, which pushes every ingroup tip
outwards by the same large constant and leaves the accessions in a thin
ring around an empty disc. Both root branches are cut back so the
ingroup fills the plot; distances *within* the ingroup are untouched,
and the caption says the branch is shortened. The rectangular figure
does not do this — it is the one that keeps every length honest.

## Reading the result

The root sits outside the Caucasian accessions and orients the tree from
the deepest ingroup split towards more recent ones. Bootstrap is how
often a clade survives resampling of the sites, not the probability that
it is true: 70 and above is worth interpreting, 90 and above is strong,
below 70 should be discussed as uncertain. Where a K group forms one
clade, the phylogeny supports the ADMIXTURE structure; mixed colours
within a clade can mean admixture, recent gene flow, close kinship — or
simply that a bifurcating tree cannot represent a reticulate history.

## Running it overnight without staying logged in

`dnaml` takes six hours or so. `scripts/finish_after_dnaml.sbatch` does
everything after it — bootstrap, validation, annotation, figures — and is
submitted with a scheduler dependency, so SLURM starts it by itself:

```
sbatch --dependency=afterany:<dnaml job id> \
       pipelines/05_tree/scripts/finish_after_dnaml.sbatch
```

`afterany`, not `afterok`, deliberately: `snphylo.vcf.sh` runs an R
plotting step *after* it has written the ML tree, so it can exit
non-zero with a perfectly good tree on disk. The script checks the tree
files rather than trusting the exit code.

It calls the same `bin/` scripts as the pipeline. `nextflow run .
--ml_tree ... --phylip ...` reproduces the same results.
