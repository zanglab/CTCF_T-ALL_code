#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF3M.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF4E 
