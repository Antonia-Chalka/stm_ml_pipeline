#!/bin/sh
#$ -N panaroo_strict
#$ -l h_rt=8:00:00
#$ -l h_vmem=4G
#$ -cwd
#$ -pe sharedmem 10
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate panaroo_env

# Set up variables, get gff files and wd
outdir="/exports/eddie/scratch/s1438773/panaroo_strict_out/" #change for test 
cd "/exports/eddie/scratch/s1438773/" #change for test 

# Print for debugging purposes
echo Outdir: "$outdir"
echo Current Working Directory: "$(pwd)"
echo annotations: $(find prokka_out/ -type f -name "*.gff") 
echo command: panaroo -i $(find prokka_out/ -type f -name "*.gff") -o $outdir -t 10 --clean-mode strict 

mkdir -p $outdir

#TODO-CHANGE COMMAND
panaroo -i $(find prokka_out/ -type f -name "*.gff") -o $outdir -t 10 --clean-mode strict

# https://gtonkinhill.github.io/panaroo/#/gettingstarted/params
# -i input files. must be gff, can also be file of gff files
# - o output folder
# --mode strict/moderate/sensitive
# --merge_paralogs


