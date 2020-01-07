#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATG10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG101 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATG16L1 
