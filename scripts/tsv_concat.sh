#!/bin/sh

categories=( bovine bovineenv bovineproducts human poultry poultryenv poultryproducts swine swineproducts)
for category in ${categories[@]}
do
    awk 'FNR>1' $(find /exports/eddie/scratch/s1438773/quast_results/$category -type f -name "transposed_report.tsv") >quast_$category.tsv
done

# Bovine all
cat quast_bovine.tsv quast_bovineenv.tsv quast_bovineproducts.tsv >quast_bovine_all.tsv

# Poultry all
cat quast_poultry.tsv quast_poultryenv.tsv quast_poultryproducts.tsv >quast_poultry_all.tsv

# Swine all
cat quast_swine.tsv quast_swineproducts.tsv >quast_swine_all.tsv

#Everything
cat *_all.tsv >quast_all.tsv

# Have to manually add human to quast_all + headers for each document
