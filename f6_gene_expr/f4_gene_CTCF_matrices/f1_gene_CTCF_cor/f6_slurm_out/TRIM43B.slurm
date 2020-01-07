#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRIM43B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM43B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIM48 
