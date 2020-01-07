#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIGD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIGIT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIMD4 
