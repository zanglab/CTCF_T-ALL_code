#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL6R.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL6R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL6ST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL7R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL9R 
