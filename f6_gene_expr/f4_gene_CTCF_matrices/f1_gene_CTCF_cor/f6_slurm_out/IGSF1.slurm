#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGSF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGSF23 
