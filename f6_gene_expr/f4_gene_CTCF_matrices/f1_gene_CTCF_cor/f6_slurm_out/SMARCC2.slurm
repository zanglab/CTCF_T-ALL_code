#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMARCC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMARCE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC1A 
