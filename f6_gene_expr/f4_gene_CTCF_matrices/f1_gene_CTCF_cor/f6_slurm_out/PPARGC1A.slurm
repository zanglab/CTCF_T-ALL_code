#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPARGC1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPARGC1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPARGC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPBPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPCDC 
