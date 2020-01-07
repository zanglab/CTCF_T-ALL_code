#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VMA21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VMA21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VMAC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VMO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VN1R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VN1R10P 
