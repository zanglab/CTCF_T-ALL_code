#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DISP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DISP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIXDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DKC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DKK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DKK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DKK3 
