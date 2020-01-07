#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VSTM2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VSTM2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSTM2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSTM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSTM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSX2 
