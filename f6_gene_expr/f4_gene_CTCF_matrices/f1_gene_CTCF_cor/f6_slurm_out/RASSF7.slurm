#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RASSF7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAVER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAVER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAX 
