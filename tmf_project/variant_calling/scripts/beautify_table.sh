#!/bin/bash
set -euo pipefail

SNPS_TXT="$1"
INDELS_TXT="$2"
GENE="$3"
SAMPLE1="$4"
SAMPLE2="$5"

KNOWN_FMF="rs61752717 rs28940579 rs28940580 rs3743930 rs11466024 rs11466045 rs11466023 rs61732874 rs104895097 rs104895083"

make_readable() {
    local input="$1"
    local output="$2"

    local ncol
    ncol=$(head -1 "$input" | awk -F'\t' '{print NF}')
    local s1=$((ncol - 1))
    local s2=$ncol

    awk -F'\t' -v s1="$s1" -v s2="$s2" -v sn1="$SAMPLE1" -v sn2="$SAMPLE2" '
    BEGIN { OFS="\t" }
    NR==1 {
        print "Chr", "Pos", "Ref", "Alt", "Gene", "Region", "Effect", "AAChange", "rsID", sn1, sn2
        next
    }
    {
        gt1 = $s1; sub(/:.*/, "", gt1)
        gt2 = $s2; sub(/:.*/, "", gt2)
        print $1, $2, $4, $5, $7, $6, $9, $10, $11, gt1, gt2
    }' "$input" | column -t -s $'\t' > "$output"
}

make_readable "$SNPS_TXT" "readable_snps.txt"
make_readable "$INDELS_TXT" "readable_indels.txt"

{
    echo "====================================="
    echo "  Variants in gene: $GENE"
    echo "====================================="

    for label_file in "SNPs:$SNPS_TXT" "Indels:$INDELS_TXT"; do
        label="${label_file%%:*}"
        file="${label_file#*:}"

        echo ""
        echo "--- ${label} ---"
        echo ""

        ncol=$(head -1 "$file" | awk -F'\t' '{print NF}')
        s1=$((ncol - 1))
        s2=$ncol

        matches=$(grep -i "$GENE" "$file" || true)

        if [ -z "$matches" ]; then
            echo "  No ${GENE} variants found."
        else
            echo "$matches" | awk -F'\t' -v s1="$s1" -v s2="$s2" \
                -v sn1="$SAMPLE1" -v sn2="$SAMPLE2" -v known="$KNOWN_FMF" '
            BEGIN {
                n = split(known, known_arr, " ")
            }
            {
                gt1 = $s1; sub(/:.*/, "", gt1)
                gt2 = $s2; sub(/:.*/, "", gt2)

                printf "  %s:%s %s>%s\n", $1, $2, $4, $5
                printf "    Region:   %s\n", $6
                printf "    Effect:   %s\n", $9
                printf "    AAChange: %s\n", $10
                printf "    rsID:     %s\n", $11
                printf "    %-10s %s\n", sn1":", gt1
                printf "    %-10s %s\n", sn2":", gt2

                is_known = 0
                for (i = 1; i <= n; i++) {
                    if ($11 == known_arr[i]) { is_known = 1; break }
                }
                if (is_known) printf "    >>> KNOWN FMF MUTATION <<<\n"

                printf "\n"
            }'
        fi
    done
} > "${GENE}_report.txt"