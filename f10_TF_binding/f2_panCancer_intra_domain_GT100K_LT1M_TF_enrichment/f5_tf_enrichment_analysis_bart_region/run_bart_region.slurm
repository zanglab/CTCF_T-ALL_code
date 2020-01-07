#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab

#Run program



infile="../f2_high_cor_high_logFC_increased_decreased/CRC_gained_center_hicor.bed"
bart2 region -i ${infile} -s hg38 -c 4 --outdir bart2_out

infile="../f2_high_cor_high_logFC_increased_decreased/T-ALL_gained_center_hicor.bed"
bart2 region -i ${infile} -s hg38 -c 4 --outdir bart2_out

