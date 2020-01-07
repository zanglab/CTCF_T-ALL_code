#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ING5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ING5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INGX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INHBA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INHBB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INHBC 
