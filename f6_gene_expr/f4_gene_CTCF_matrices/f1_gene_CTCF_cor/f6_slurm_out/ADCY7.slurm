#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADCY7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCYAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCYAP1R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADD1 
