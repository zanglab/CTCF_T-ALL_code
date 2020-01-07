#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHRDL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRDL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRDL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRFAM7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHRM3 
