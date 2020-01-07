#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPP30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPP40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPRD1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPRD1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPRD2 
