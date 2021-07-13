#!/bin/sh
#$ -N prokka_bovine_env
#$ -cwd
#$ -l h_rt=48:00:00 
#$ -l h_vmem=5G
#$ -pe sharedmem 8 
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules (needed?)
. /etc/profile.d/modules.sh

# anaconda env
module load anaconda
source activate prokka_env

hosts=( bovineenv )

for host in "${hosts[@]}"
do
    echo Host: "$host"
    outdir=/exports/eddie/scratch/s1438773/prokka_out/$host
    mkdir -p "$outdir"
    echo Outdir: "$outdir"

    # Get all assembly files of a host into an array
    cd /exports/eddie/scratch/s1438773/assemblies/STm/usa/$host || exit
    echo Current Working Directory: $(pwd)

    assemblies=( * )
    echo assemblies: "${assemblies[@]}"

    for fullfile in "${assemblies[@]}"
    do
          assembly=$(basename "$fullfile" .fasta)
          outhostdir=$outdir/$assembly
          mkdir -p "$outhostdir"
          echo Assembly: "$assembly"
          echo Assembly file and path: "$fullfile"
          echo Output directory: "$outhostdir"
          
          echo prokka --outdir $outhostdir --prefix $assembly --force --proteins /exports/cmvm/eddie/eb/groups/gally_grp/annita/stm_proteinref.fasta --centre X --compliant "$fullfile"
          prokka --outdir $outhostdir --prefix $assembly --force  --proteins /exports/cmvm/eddie/eb/groups/gally_grp/annita/stm_proteinref.fasta --centre X --compliant "$fullfile"
    done
done
