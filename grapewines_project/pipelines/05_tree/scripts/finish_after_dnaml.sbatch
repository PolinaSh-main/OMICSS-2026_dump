#!/bin/bash
#SBATCH --job-name=tree_finish
#SBATCH --output=/path/to/workdir/tree_rooted/logs/finish_%j.out
#SBATCH --error=/path/to/workdir/tree_rooted/logs/finish_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --partition=compute

#
# Everything after dnaml: bootstrap, validation, tip annotation, figures.
#
# Submitted with --dependency=afterany on the dnaml job, so SLURM starts
# it by itself and nothing has to stay logged in. "afterany" rather than
# "afterok" on purpose: snphylo.vcf.sh runs an R plotting step *after*
# it has already written the ML tree, so the script can exit non-zero
# with a perfectly good tree on disk. The tree files are checked here
# instead of trusting the exit code.
#
# These are the same bin/ scripts that pipelines/05_tree calls. Running
# them directly avoids a Nextflow driver inside a compute-node
# allocation; `nextflow run . --ml_tree ... --phylip ...` reproduces
# this from the repository.
#

set -uo pipefail

REPO=/path/to/workdir
PROJ=$REPO/grapewines_project
BIN=$PROJ/pipelines/05_tree/bin
RUN=$REPO/tree_rooted
OUT=$PROJ/results/tree
SNP=/path/to/shared-data/phylogeny/software/snphylo

export PATH=/path/to/tools/bin:$PATH
export R_LIBS=$SNP/R_LIBS

ML=$RUN/results/cauc_rooted.ml.tree
PHY=$RUN/results/cauc_rooted.phylip.txt
BS=$OUT/snphylo/cauc_rooted.bs.tree

CORES=${SLURM_CPUS_PER_TASK:-8}

step () { echo; echo "=== $* :: $(date -Is) ==="; }

fail () { echo; echo "FAILED: $*"; exit 1; }


step "checking what dnaml left"

for f in "$ML" "$PHY"; do
    [ -s "$f" ] || fail "missing or empty: $f -- read tree_rooted/logs/ml_*.err"
    echo "  ok  $f  ($(stat -c%s "$f") bytes)"
done

echo "  tips in ML tree: $(tr ',' '\n' < "$ML" | grep -c ':')"
echo "  sites          : $(head -1 "$PHY" | awk '{print $2}')"


mkdir -p "$OUT/snphylo" "$OUT/plots" || fail "cannot create $OUT"

for f in cauc_rooted.ml.tree cauc_rooted.ml.txt cauc_rooted.phylip.txt \
         cauc_rooted.fasta cauc_rooted.id.txt; do
    [ -s "$RUN/results/$f" ] && cp -f "$RUN/results/$f" "$OUT/snphylo/"
done


step "bootstrap, 100 replicates on $CORES cores"

"$BIN/bootstrap_tree.R" \
    --tree "$ML" \
    --phylip "$PHY" \
    --replicates 100 \
    --cores "$CORES" \
    --out "$BS" || fail "bootstrap"


step "validation"

"$BIN/validate_tree.py" \
    --ml-tree "$ML" \
    --bs-tree "$BS" \
    --expect-tips 413 \
    --outgroup ZZ01 \
    --out "$OUT/tree_validation.txt" || fail "validation"


step "tip annotation from the K=7 Q matrix"

"$BIN/annotate_tips.py" \
    --tree "$BS" \
    --q "$PROJ/results/admixture/Q/cauc_filtered.final.7.Q" \
    --fam "$PROJ/reference/cauc_filtered.final.fam" \
    --metadata "$PROJ/reference/cauc_grape_metadata.csv" \
    --min-q 0.75 \
    --outgroup ZZ01 \
    --outdir "$OUT" || fail "annotation"


step "figures"

for layout in fan rectangular; do

    "$BIN/plot_tree.py" \
        --tree "$BS" \
        --annotation "$OUT/tip_annotation.tsv" \
        --labels "$OUT/group_labels.tsv" \
        --layout "$layout" \
        --outgroup ZZ01 \
        --min-support 70 \
        --outdir "$OUT/plots" || fail "$layout figure"

done


step "run summary"

{
    echo "input           /path/to/shared-data/phylogeny/rooted/input/cauc_rooted_tree_input.413samples.vcf.gz"
    echo "samples         413 (412 Caucasian + ZZ01)"
    echo "variants in     119195"
    echo "sites retained  $(head -1 "$PHY" | awk '{print $2}')"
    echo "maf             >= 0.10"
    echo "missing rate    <= 0.10"
    echo "ld threshold    0.10"
    echo "chromosomes     1-19"
    echo "low-depth screen skipped (-r)"
    echo "outgroup        ZZ01"
    echo "method          maximum likelihood, PHYLIP dnaml via SNPhylo"
    echo "snphylo         $SNP/snphylo.vcf.sh"
    echo "bootstrap       100 replicates, phangorn::bootstrap.pml, optNni"
    echo "slurm jobs      dnaml ${DNAML_JOB:-unknown}, rest $SLURM_JOB_ID"
    echo
    echo "ml tree         $OUT/snphylo/cauc_rooted.ml.tree"
    echo "bootstrap tree  $BS   <- the final tree"
    echo "annotation      $OUT/tip_annotation.tsv"
    echo "figures         $OUT/plots/"
} > "$OUT/run_summary.txt"

cat "$OUT/run_summary.txt"


step "done"

ls -l "$OUT" "$OUT/snphylo" "$OUT/plots"

echo
echo "Everything is under $OUT"
echo "Commit it with:  cd $REPO && git add grapewines_project/results/tree && git commit"
