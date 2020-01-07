#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OVOL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OVOL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OVOL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXA1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXCT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXCT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXER1 
