#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SP9 
