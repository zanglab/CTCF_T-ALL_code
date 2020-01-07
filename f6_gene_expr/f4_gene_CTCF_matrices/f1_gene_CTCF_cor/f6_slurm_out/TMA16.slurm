#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMA16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMA16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMBIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMBIM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMBIM6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC1 
