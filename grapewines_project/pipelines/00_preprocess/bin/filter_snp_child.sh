#SCRIPT FOR SEPERATION OF SNPs AND INDELS
# NOTE: Run this script from the directory where the "log" directory is located,
#       Example: /mnt/nas1/proj/omicss26/admixture
#
# PURPOSE:
#   This script performs snp filtering.
#
# PARAMETERS:
#   1: input_dir - Directory where vcf/gvcf files are located.
#   2: store_plink - Directory where plink binary format files will be stored (we would have to convert vcf into plink binary files)
#   3: output_dir - Directory where filtered result will be located
# SAMPLE USAGE:
#   In a parent script: src/filter_snp/filter_snp_parent.sh <input_dir> > <output_dir>
#
# IMPORTANT:
#   - Run from a parent script.

# Check if the correct number of arguments is provided

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file> <store_plink_dir> <output_dir>"
    exit 1
fi

input_file="$1"
store_plink="$2"
output_dir="$3"
plink="/mnt/nas0/proj/vine/user_projects/shengchang/soft/plink/plink"

mkdir -p "$store_plink"
mkdir -p "$output_dir"

set -e

# Step1: Filtering part until LD prunning

echo "Doing filtering (no LD prunning)"

${plink} --vcf "${input_file}" \
         --make-bed \
         --double-id \
         --keep-allele-order \
         --allow-extra-chr 0 \
         --chr-set 19 \
         --biallelic-only strict \
         --maf 0.05 \
         --geno 0 \
        --out "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0"


echo "Done first part of filtering"

# sort the bim file for ld prunning

echo "Sorting bim file for LD prunning"
sort -k1,1n -k4,4n "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.bim" > "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.bim.sorted"
echo "Done sorting"

cut -f2 "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.bim.sorted" > variant_order.txt

echo "Order of variants is stored"

echo "Sorting  plink binary formal files"
${plink} --bfile "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0" \
      --extract variant_order.txt \
      --make-bed \
      --out "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.sorted"
echo "Sorting completely done!"

#STEP 6: LD pruning

echo "Doing LD prunning"
${plink} --bfile "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.sorted" \
         --indep-pairwise 50 5 0.4 \
         --out "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.ld_pruned"

echo "LD pruning done (window=50, step=5, r�=0.04)."


echo "Filtering complete!"

echo "Making final filtered file"

${plink} --bfile "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.sorted" \
         --extract "${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.ld_pruned.prune.in" \
         --make-bed \
         --out "${output_dir}/cauc_filtered.final"

echo "Filtering compleletly completed!"
  
  
echo "filtering vcf out of vcf"


# THIS PARRT TAAKES FILTERED VARIANTS FROM ORIGINAL VCF AND GIVES FILTERED VCF 

#bcftools view -Oz -o "${output_dir}/cauc_filtered.final.vcf.gz" -i "ID=@${store_plink}/cauc_plink.1-19.maf0.05.bi.g0.ld_pruned.prune.in" "${store_plink}/cauc_grape.filtered.subset.snp.vcf.gz"
# Index the filtered VCF
#tabix -p vcf "${output_dir}/cauc_filtered.final.vcf.gz"

#echo "Done!"






































