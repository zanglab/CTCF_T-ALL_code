#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM38 
