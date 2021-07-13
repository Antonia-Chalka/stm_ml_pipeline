#!/bin/sh
#$ -N quast_swine
#$ -cwd
#$ -l h_rt=10:00:00 
#$ -l h_vmem=2G
#$ -pe scatter 5

# Initialise the environment modules (needed?)
. /etc/profile.d/modules.sh

# anaconda env
module load anaconda
source activate quast


hosts=( swine swineproducts )
quast_dir=/exports/cmvm/eddie/eb/groups/gally_grp/annita/software/quast-5.0.2

for host in "${hosts[@]}"
do
    echo Host: "$host"
    outdir=/exports/eddie/scratch/s1438773/quast_results/$host
    mkdir -p "$outdir"
    echo Outdir: "$outdir"

    # Get all assembly files of a host into an array
    cd "/exports/cmvm/eddie/eb/groups/gally_grp/annita/assemblies/STm/usa/$host" || exit
    echo Current Working Directory: `pwd`

    assemblies=( * )
    echo assemblies: "${assemblies[@]}"

    #output array to file to document which assemblies got inputted
    for fullfile in "${assemblies[@]}"
    do
          assembly=$(basename "$fullfile" .fasta)
          outhostdir=$outdir/$assembly
          mkdir -p "$outhostdir"
          echo Assembly: "$assembly"
          echo Assembly file and path: "$fullfile"
          echo Output directory: "$outhostdir"
          
          echo Quast command: python $quast_dir/quast.py -o "$outhostdir" -r /exports/cmvm/eddie/eb/groups/gally_grp/annita/assemblies/STm/NC_003198.1_quast_ref.fasta  "$fullfile"
          python $quast_dir/quast.py -o "$outhostdir" -r /exports/cmvm/eddie/eb/groups/gally_grp/annita/assemblies/STm/NC_003198.1_quast_ref.fasta  "$fullfile"

    done
done
