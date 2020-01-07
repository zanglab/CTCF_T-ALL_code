#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRPAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRPAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRPPRC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC10B 
