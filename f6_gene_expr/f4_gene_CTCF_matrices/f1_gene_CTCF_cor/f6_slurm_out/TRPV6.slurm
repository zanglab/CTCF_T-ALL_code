#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRPV6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPV6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRRAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRUB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSACC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSC1 
