#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZSCAN31.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZSCAN9 
