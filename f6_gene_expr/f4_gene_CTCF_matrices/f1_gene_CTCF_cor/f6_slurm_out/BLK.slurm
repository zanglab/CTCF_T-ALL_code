#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BLK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BLK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLMH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLNK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S2 
