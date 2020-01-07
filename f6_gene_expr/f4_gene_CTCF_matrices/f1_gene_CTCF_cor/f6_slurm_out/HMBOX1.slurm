#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HMBOX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HMBOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMBS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMCES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMG20A 
