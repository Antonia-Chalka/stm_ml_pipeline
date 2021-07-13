#!/bin/sh
#$ -N snippy_bov_env
#$ -cwd
#$ -l h_rt=48:00:00 
#$ -l h_vmem=3G
#$ -pe sharedmem 10
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
module load java
source activate snippy_env

ref_seq="/exports/eddie/scratch/s1438773/panaroo_moderate_out/pan_genome_reference.fa"
hosts=( bovineenv )

for host in "${hosts[@]}"
do
    echo Host: "$host"
    
    #change for test
    outdir=/exports/eddie/scratch/s1438773/snippy_out/$host
    mkdir -p "$outdir"
    echo Outdir: "$outdir"

    # Get all assembly files of a host into an array 
    #change for test
    cd /exports/eddie/scratch/s1438773/assemblies_qc/STm/usa/$host || exit
    echo Current Working Directory: $(pwd)
    assemblies=( * )
    echo assemblies: "${assemblies[@]}"

    #output array to file to document which annotations got inputted
    for fullfile in "${assemblies[@]}"
    do
          assembly=$(basename "$fullfile" .fasta)
          outhostdir=$outdir/$assembly
          mkdir -p "$outhostdir"
          echo Assembly: "$assembly"
          echo Assembly file and path: "$fullfile"
          echo Output directory: "$outhostdir"
          
          echo snippy command: snippy --report --cpus 10 --ram 29 --outdir  "$outhostdir" --prefix "$assembly" --ref  "$ref_seq" --ctgs "$fullfile" --force
          snippy command: snippy --report --cpus 10 --ram 29 --outdir  "$outhostdir" --prefix "$assembly" --ref  "$ref_seq" --ctgs "$fullfile" --force
          
    done
done

# get folders folders=`find . -mindepth 2 -type d`
# echo ${folders[@]}
# snippy-core --ref refseq ${folders[@]}
