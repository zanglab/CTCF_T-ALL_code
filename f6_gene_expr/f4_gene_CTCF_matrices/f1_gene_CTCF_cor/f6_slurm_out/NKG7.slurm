#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NKG7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NKG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKIRAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKIRAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKTR 
