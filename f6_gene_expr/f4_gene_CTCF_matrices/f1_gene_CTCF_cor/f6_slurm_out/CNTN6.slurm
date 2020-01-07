#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNTN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP3P2 
