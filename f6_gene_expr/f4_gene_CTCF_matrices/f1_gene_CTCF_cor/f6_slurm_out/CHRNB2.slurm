#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHRNB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRND 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRNG 
