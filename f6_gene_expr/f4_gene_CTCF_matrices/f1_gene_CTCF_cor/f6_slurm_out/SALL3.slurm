#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SALL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SALL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SALL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD12 
