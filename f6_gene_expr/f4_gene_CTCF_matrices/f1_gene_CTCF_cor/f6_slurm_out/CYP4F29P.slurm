#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP4F29P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F29P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F30P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F35P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F62P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4F8 
