#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDKN1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN2AIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN2AIPNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDKN2C 
