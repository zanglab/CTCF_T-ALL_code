#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CFAP20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP221 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CFAP45 
