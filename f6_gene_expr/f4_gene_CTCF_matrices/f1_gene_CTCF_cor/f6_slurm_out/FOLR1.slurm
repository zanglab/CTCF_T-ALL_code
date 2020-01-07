#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOLR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOLR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOLR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOLR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOPNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOSB 
