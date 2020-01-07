#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STARD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STARD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAT1 
