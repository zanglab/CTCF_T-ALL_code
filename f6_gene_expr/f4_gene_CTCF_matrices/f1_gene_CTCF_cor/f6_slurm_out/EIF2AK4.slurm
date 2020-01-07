#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF2AK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2AK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2B5 
