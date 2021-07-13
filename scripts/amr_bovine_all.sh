#!/bin/sh
#$ -N amr_bovine_all
#$ -cwd
#$ -l h_rt=40:00:00
#$ -l h_vmem=650M
#$ -pe sharedmem 3
#$ -m beas
#$ -M s1438773@ed.ac.uk
#$ -hold_jid prokka_bovine_an

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate amrfinder_env

hosts=( bovine bovineproducts bovineenv)

for host in "${hosts[@]}"
do
    outdir=/exports/eddie/scratch/s1438773/amr_out/$host # Change for test
    mkdir -p "$outdir"
    cd "/exports/eddie/scratch/s1438773/prokka_out/$host" || exit # Input dir - Change for test
    annotations=$(find * -type f -name "*.faa") # Get all annotation files of a host into an array
    
    echo Host: "$host"
    echo Outdir: "$outdir"
    echo Current Working Directory: "$(pwd)"
    echo annotations: "$annotations"
    
    for annotation in $annotations
    do
        filename=$(basename "$annotation" .faa)
        
        echo annotation: "$filename" path: "$annotation" 
        echo Output directory: "$outdir"
        echo amrfinder command: amrfinder --plus -O Salmonella -p "$annotation" -o "$outdir"/"$filename"_amr.tsv
        
        amrfinder --plus -O Salmonella -p "$annotation" -o "$outdir"/"$filename"_amr.tsv
    done
done
