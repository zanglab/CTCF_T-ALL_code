#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CACNB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG3 
