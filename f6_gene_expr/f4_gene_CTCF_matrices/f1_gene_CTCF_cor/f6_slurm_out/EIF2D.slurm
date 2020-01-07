#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EIF2D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2S2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF2S3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EIF3B 
