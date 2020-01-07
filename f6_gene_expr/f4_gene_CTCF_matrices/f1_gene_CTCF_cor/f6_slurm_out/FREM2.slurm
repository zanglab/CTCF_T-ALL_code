#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FREM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FREM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FREM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG1BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG1CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG1DP 
