#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELOVL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL7 
