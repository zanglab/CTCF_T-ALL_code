#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BLOC1S3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S5-TXNDC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLOC1S6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLVRA 
