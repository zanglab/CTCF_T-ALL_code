#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPSG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPTE2 
