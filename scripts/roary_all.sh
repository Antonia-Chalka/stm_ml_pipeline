#!/bin/sh
#$ -N roary_all
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=6G
#$ -pe sharedmem 10
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate roary_env

# Set up variables, get gff files and wd
outdir="/exports/eddie/scratch/s1438773/roary_out_nuc" #change for test
cd "/exports/eddie/scratch/s1438773/prokka_out/" #change for test

# Print for debugging purposes
echo Outdir: "$outdir"
echo Current Working Directory: "$(pwd)"
echo annotations: $(find -type f -name "*.gff")
echo roary command: roary -p 10 -v -s -f $outdir -o "stm_all_clustered_proteins" -r $(find -type f -name "*.gff")

roary -a # Check Roary dependencies
roary -p 10 -v -s -f $outdir -o "stm_all_clustered_proteins" -r $(find -type f -name "*.gff")

#Usage:   roary [options] *.gff
#         -a        check dependancies and print versions
#         -p INT    number of threads [1]
#         -v        verbose output to STDOUT
#         -s        dont split paralogs
#         -f STR    output directory [.]
#         -o STR    clusters output filename [clustered_proteins]
#         -qc       generate QC report with Kraken
#         -r        create R plots, requires R and ggplot2

#         -e        create a multiFASTA alignment of core genes using PRANK
#         -n        fast core gene alignment with MAFFT, use with -e (faster but less accurate, optional)
#         -i        minimum percentage identity for blastp [95]
#         -cd FLOAT percentage of isolates a gene must be in to be core [99]
#         -ap       allow paralogs in core alignment
#         -z        dont delete intermediate files
#         -y        add gene inference information to spreadsheet, doesnt work with -e
#         -iv STR   Change the MCL inflation value [1.5]
# As I'm using this only for piggy input, dont need to generate pan genome
