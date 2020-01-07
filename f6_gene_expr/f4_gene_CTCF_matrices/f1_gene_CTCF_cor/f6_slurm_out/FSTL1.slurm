#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FSTL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FSTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSTL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSTL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSTL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTCDNL1 
