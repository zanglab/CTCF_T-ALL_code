#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PARPBP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PARPBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARVA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARVB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PARVG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PASD1 
