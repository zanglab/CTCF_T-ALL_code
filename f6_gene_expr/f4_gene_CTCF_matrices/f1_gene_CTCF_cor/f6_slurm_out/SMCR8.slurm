#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMCR8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMCR8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMDT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMG1P3 
