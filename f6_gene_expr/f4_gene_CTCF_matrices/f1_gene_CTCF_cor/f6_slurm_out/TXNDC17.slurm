#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TXNDC17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNIP 
