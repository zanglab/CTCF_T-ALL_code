#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FATE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAU 
