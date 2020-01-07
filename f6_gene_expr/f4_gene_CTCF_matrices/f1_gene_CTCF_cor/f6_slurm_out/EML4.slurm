#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EML4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EML4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EML5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EML6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMP3 
