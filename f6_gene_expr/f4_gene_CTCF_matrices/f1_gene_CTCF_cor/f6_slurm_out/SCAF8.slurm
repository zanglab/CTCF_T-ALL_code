#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCAF8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAMP4 
