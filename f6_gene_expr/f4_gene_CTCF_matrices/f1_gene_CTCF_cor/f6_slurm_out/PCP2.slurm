#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCP4L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK1N 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK2 
