#!/bin/sh
#$ -N panaroo_msa_core
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=3G 
#$ -pe sharedmem 32
#$ -m beas
#$ -M s1438773@ed.ac.uk

# Initialise the environment modules 
. /etc/profile.d/modules.sh
module load anaconda
source activate panaroo_env 

outdir="/exports/eddie/scratch/s1438773/panaroo_moderate_out/"
cd $outdir

echo Current Working Directory: "$(pwd)"
echo  command: panaroo-msa -t 32 -a core -o $outdir

panaroo-msa --version
panaroo-msa -t 32 -a core -o $outdir

#optional arguments:
#  -h, --help            show this help message and exit
#  -t N_CPU, --threads N_CPU
#                        number of threads to use (default=1)
#  --verbose             print additional output
#  --version             show program's version number and exit
#Input/output:
#  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
#                        location of the Panaroo output directory
#Gene alignment:
#  -a {core,pan}, --alignment {core,pan}
#                        Output alignments of core genes or all genes. Options
#                        are 'core' and 'pan'. Default: 'None'
#  --aligner {clustal,prank,mafft}
#                        Specify an aligner. Options:'prank', 'clustal', and
#                        default: 'mafft'
#  --core_threshold CORE
#                        Core-genome sample threshold (default=0.95)
