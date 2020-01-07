#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TXN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TXN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXNDC16 
