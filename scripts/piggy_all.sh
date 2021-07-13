#!/bin/sh
#$ -N piggy_all
#$ -cwd
#$ -l h_rt=40:00:00
#$ -l h_vmem=3G
#$ -pe sharedmem 6
#$ -m beas
#$ -M s1438773@ed.ac.uk
#$ -hold_jid roary_all_nuc

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate piggy_env

# Set up variables
outdir="/exports/eddie/scratch/s1438773/piggy_out/" 
roarydir="/exports/eddie/scratch/s1438773/roary_out_nuc/" 
inputdir="/exports/eddie/scratch/s1438773/piggy_in/" #TODO-MAKEDIR
piggydir="/exports/cmvm/eddie/eb/groups/gally_grp/annita/software/piggy/bin"

# Print for debugging purposes
echo Current Working Directory: "$(pwd)"
echo  command: "$piggydir"/piggy -R on -t 6 -r $roarydir -o $outdir -i $inputdir

mkdir -p $outdir
"$piggydir"/piggy -R on -t 6 -r $roarydir -o $outdir -i $inputdir

#--in_dir|-i  <STR>  input folder [default - current folder]
#--out_dir|-o  <STR>  output folder [default - current folder/piggy_out]
#--roary_dir|-r  <STR>  folder where roary output is stored [required]
#--threads|-t <INT> threads [default - 1]
#--nuc_id|-n <INT> min percentage nucleotide identity [default - 90]
#--len_id|-l <INT> min percentage length identity [default - 90]
#--edges|-e  keep IGRs at the edge of contigs [default - off]
#--size|-s <STR> size of IGRs to extract [i-j] [default 30-1000]
#--method|-m <STR> method for detecting switched IGRs [g - gene_pair, u - upstream] [default - g]
#--R_plots|-R  make R plots (requires R, Rscript, ggplot2, reshape2) [default - off]
#--fast|-f  fast mode (doesn't align IGRs or detect switched regions) [default - off]
#--help|-h  help
#--version|-v  version

