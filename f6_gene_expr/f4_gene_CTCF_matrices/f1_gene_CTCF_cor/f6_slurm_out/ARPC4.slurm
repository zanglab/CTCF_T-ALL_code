#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARPC4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC4-TTLL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPC5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARPP19 
