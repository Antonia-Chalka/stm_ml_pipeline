#!/bin/sh
#$ -N scoary_collapsed
#$ -cwd
#$ -l h_rt=3:00:00 
#$ -l h_vmem=8G 
#$ -pe sharedmem 6 
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate scoary_env

# Set up variables, get gff files and wd
traitfile="/exports/cmvm/eddie/eb/groups/gally_grp/annita/scoary_hostdata.csv"
genefile="/exports/eddie/scratch/s1438773/panaroo_moderate_out/gene_presence_absence_roary.csv"

outdir="/exports/eddie/scratch/s1438773/scoary_out/collapsed" #change for test #TODO-CHANGE
mkdir -p $outdir

# Print for debugging purposes
echo  command: scoary -t $traitfile -g$genefile -o $outdir --no-time --threads 5 --collapse -p 1.0

#TODO-CHANGE COMMAND
scoary --version
scoary -t $traitfile -g $genefile -o $outdir --no-time --threads 5 --collapse -p 1.0


# --collapse flag
# Adding this flag to the command line will collapse genes that are identically distributed in your sample. 
# For example, plasmid genes are often inherited together and as such will not add any information individually. 
#From a statistical point of view, this is more correct, as in the opposite case the program will test (and correct for) multiple identical null hypotheses, 
# thus unfairly penalizing your results by multiple comparisons correction.
# may use it as would do sth like that either way in building the RF model?
# but if i report for every gene then dont really need it...
# run scoary 2 times, one collapsed, one not?

# can set p cuttooff to 1 to get every gene, might be seful as i can see the distribution of p values? plus i can easily filter them myself in r

