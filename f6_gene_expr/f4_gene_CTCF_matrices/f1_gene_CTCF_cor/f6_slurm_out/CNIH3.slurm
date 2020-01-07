#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNIH3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNIH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNIH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNKSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNKSR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNKSR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNN1 
