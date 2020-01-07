#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSSK6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSSK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSTA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSTD3 
