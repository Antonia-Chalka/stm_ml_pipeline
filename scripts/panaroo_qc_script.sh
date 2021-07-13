#!/bin/sh
#$ -N panaroo_qc 
#$ -cwd
#$ -l h_rt=40:00:00
#$ -l h_vmem=10G
#$ -pe sharedmem 5
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate panaroo_env

# Set up variables, get gff files and wd
outdir="/exports/eddie/scratch/s1438773/panaroo_qc_out/" #change for test 
cd "/exports/eddie/scratch/s1438773/" #change for test

# Print for debugging purposes
echo Outdir: "$outdir"
echo Current Working Directory: "$(pwd)"
echo annotations: $(find prokka_out/ -type f -name "*.gff") 
echo command: panaroo-qc -i $(find prokka_out/ -type f -name "*.gff") -o $outdir -t 5 --graph_type all --ref_db refseq.genomes.k21s1000.msh

mkdir -p $outdir
panaroo-qc -i $(find prokka_out/ -type f -name "*.gff") -o $outdir -t 5 --graph_type all --ref_db refseq.genomes.k21s1000.msh

# -i input files. must be gff, can also be file of gff files
# - o output folder
# --mode strict/moderate/sensitive
# --merge_paralogs

