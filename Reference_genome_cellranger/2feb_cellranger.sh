#!/bin/sh
#SBATCH --job-name=Cellranger
#SBATCH --mail-user=andrew.ramos@imm.ox.ac.uk
#SBATCH --mail-type=END,FAIL,BEGIN,ALL 
#SBATCH --output=%j_%x.log
#SBATCH --error=%j_%x.log
#SBATCH --ntasks=20
#SBATCH --mem=256gb

module load cellranger/6.1.2
                                cellranger count --id=D52N_miSFITs \
                                 --fastqs=/t1-data/project/tsslab/andramos/sc_rnaseq/data/scrnaseq_raw_fastq \
                                 --sample=ASR_Novo_10X-SCI7T045-AK40147_HMYHGDSX2,ASR_Novo_10X-SCI7T045-SCI5T045_HMYHGDSX2 \
                                 --transcriptome=/home/a/andramos/t1data/ASRGenome/0_references/feb2/GRCh38_D52NCAR/D52N_miSFITs \
