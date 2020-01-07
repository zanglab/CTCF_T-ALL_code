#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INSIG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INSIG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSM1 
