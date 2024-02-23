#!/bin/sh
#SBATCH --job-name=cellranger_makereference
#SBATCH --mail-user=andrew.ramos@imm.ox.ac.uk
#SBATCH --mail-type=END,FAIL,BEGIN,ALL 
#SBATCH --output=%j_%x.log
#SBATCH --error=%j_%x.log
#SBATCH --ntasks=20
#SBATCH --mem=128gb

module load cellranger/6.1.2
                                cellranger mkref --genome=D52N_miSFITs \
                                 --fasta=/home/a/andramos/t1data/ASRGenome/0_references/feb2/0_building/input_for_mkref/genome+D52NCAR.fa \
                                 --genes=/home/a/andramos/t1data/ASRGenome/0_references/feb2/0_building/input_for_mkref/genes-filtered+D52NCAR.gtf \
