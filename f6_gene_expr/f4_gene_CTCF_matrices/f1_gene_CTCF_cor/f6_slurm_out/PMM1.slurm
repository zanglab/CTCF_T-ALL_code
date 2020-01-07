#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PMM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PMM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMP22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMPCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMPCB 
