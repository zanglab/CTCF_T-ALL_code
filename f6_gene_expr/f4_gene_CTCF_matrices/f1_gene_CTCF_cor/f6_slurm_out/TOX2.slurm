#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53AIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53BP1 
