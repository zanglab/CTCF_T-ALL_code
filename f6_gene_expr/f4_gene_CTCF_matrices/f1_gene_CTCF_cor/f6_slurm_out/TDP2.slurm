#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TDP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TDP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD3 
