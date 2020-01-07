#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF5A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF5A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF5AL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELAC1 
