#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBRG1 
