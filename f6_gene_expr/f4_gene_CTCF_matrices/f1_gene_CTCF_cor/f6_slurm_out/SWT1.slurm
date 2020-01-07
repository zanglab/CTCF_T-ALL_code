#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SWT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SWT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYBU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCE1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYCE2 
