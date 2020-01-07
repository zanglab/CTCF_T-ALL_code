#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAV3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAV3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBEA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBEAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBEAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBEAP1 
