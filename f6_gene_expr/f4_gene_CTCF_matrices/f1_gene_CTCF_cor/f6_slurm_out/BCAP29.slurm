#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCAP29.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAP29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAP31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAS2 
