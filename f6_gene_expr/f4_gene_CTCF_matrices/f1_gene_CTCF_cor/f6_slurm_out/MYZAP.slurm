#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYZAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYZAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MZB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MZT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MZT2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MZT2B 
