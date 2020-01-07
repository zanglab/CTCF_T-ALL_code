#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFAP58.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP69 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP97 
