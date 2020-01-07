#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMIM8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMIM9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMKR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMLR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMN2 
