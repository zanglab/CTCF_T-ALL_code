#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSPH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSTPIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSTPIP2 
