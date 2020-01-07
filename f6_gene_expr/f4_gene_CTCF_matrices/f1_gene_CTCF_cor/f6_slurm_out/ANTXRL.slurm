#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANTXRL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANTXRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANTXRLP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANXA13 
