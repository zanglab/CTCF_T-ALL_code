#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ESR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ESR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESRRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESRRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ESRRG 
