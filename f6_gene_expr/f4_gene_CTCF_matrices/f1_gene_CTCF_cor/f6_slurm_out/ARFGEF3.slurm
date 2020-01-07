#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARFGEF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFGEF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARFRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARG2 
