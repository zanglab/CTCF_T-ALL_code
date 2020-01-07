#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTMR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR2 
