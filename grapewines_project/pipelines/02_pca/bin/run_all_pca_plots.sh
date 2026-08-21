#!/bin/bash
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH -o /path/to/workdir/grapewines_project/admixture/log/pca_all.out
#SBATCH -e /path/to/workdir/grapewines_project/admixture/log/pca_all.err


PROJECT="/path/to/workdir/grapewines_project/admixture"

EIGENVEC="$PROJECT/pca_results/cauc_pca.eigenvec"
EIGENVAL="$PROJECT/pca_results/cauc_pca.eigenval"

METADATA="/path/to/shared-data/data/metadata/cauc_grape_metadata.csv"

OUTBASE="$PROJECT/pca_all_plots"


mkdir -p "$OUTBASE"
mkdir -p "$PROJECT/log"


# Получаем имена всех колонок metadata
python3 - <<EOF > columns.txt
import pandas as pd

df = pd.read_csv("$METADATA")

for c in df.columns:
    if c != "Column 1":
        print(c)
EOF


while IFS= read -r col
do

    safe_col=$(echo "$col" | sed 's/[^A-Za-z0-9_]/_/g')

    echo "Processing: $col"

    mkdir -p "$OUTBASE/$safe_col"

    python3 "/path/to/workdir/grapewines_project/admixture/scripts/pca_plot.py" \
        --eigenvec "$EIGENVEC" \
        --eigenval "$EIGENVAL" \
        --metadata "$METADATA" \
        --color-column "$col" \
        --outdir "$OUTBASE/$safe_col"

done < columns.txt


echo "DONE"