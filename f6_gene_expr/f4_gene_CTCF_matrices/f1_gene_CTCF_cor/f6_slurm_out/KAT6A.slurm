#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KAT6A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KATNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KATNAL1 
