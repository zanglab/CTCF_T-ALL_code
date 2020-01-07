#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF3H.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3I 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3L 
