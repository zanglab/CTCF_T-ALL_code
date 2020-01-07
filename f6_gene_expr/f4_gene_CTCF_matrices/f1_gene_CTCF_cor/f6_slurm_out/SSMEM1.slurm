#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SSMEM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SSMEM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSR2 
