#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=30000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o download.out

# 
# prefetch -X 30G -v SRR2094193	#GSM1816076 H3K27ac-A5
# prefetch -X 30G -v SRR2094186 #GSM1816069 H3K4me1-A5
fastq-dump --gzip --split-3 /nv/vol190/zanglab/zw5j/data/ncbi/public/sra/SRR2094193.sra -O ./
fastq-dump --gzip --split-3 /nv/vol190/zanglab/zw5j/data/ncbi/public/sra/SRR2094186.sra -O ./
mv SRR2094193.fastq.gz erythroid_progenitor_cells_GSM1816076_H3K27ac.fastq.gz
mv SRR2094186.fastq.gz erythroid_progenitor_cells_GSM1816069_H3K4me1.fastq.gz


# prefetch -X 30G -v SRR6987861 #GSM3094374 EOL1_H3K27Ac
# prefetch -X 30G -v SRR6987851 #GSM3094366 EOL1_input
fastq-dump --gzip --split-3 /nv/vol190/zanglab/zw5j/data/ncbi/public/sra/SRR6987861.sra -O ./
fastq-dump --gzip --split-3 /nv/vol190/zanglab/zw5j/data/ncbi/public/sra/SRR6987851.sra -O ./ &
mv SRR6987861.fastq.gz EOL1_H3K27Ac.fastq.gz
mv SRR6987851.fastq.gz EOL1_input.fastq.gz


