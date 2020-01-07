#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTAFR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTAFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCD1 
