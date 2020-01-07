#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDHX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDHX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDIA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDIA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDIA3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDIA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDIA5 
