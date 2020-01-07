#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CER1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CERCAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CERK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CERKL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CERS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CERS2 
