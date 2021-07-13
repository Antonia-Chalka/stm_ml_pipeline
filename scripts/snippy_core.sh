#!/bin/sh
#$ -N snippy_core
#$ -cwd
#$ -l h_rt=5:00:00 
#$ -l h_vmem=4G
#$ -pe sharedmem 20
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load java
module load anaconda
source activate snippy_env

# Set up directory and file paths
ref_seq="/exports/eddie/scratch/s1438773/panaroo_moderate_out/pan_genome_reference.fa"
snippy_dir="/exports/eddie/scratch/s1438773/snippy_out"

# Make and switch to output directory
outdir="/exports/eddie/scratch/s1438773/snippy_core_out/"
mkdir -p "$outdir"
cd "$outdir" || exit

# Get all snippy folders from snippy_out
folders=`find ${snippy_dir} -maxdepth 2 -mindepth 2 -type d`
echo ${folders[@]}

# Checks
snippy-core --version
snippy-core --check

# Snippy command
echo snippy command: snippy-core --ref  "$ref_seq" ${folders[@]}
snippy-core --ref  "$ref_seq" ${folders[@]}

# will i need to use ram/cpu check? tested snippy core with --cpus and it does not seem to have that option so uh????
# --cpus 30 --ram 60

#source activate fasttree_env
# for fasttree use core.aln
#fasttree -out core_tree.newick -nt core.aln 


