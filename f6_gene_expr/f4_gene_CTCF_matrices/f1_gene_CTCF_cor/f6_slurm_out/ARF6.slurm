#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARF6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFGAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFGAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFGEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFGEF2 
