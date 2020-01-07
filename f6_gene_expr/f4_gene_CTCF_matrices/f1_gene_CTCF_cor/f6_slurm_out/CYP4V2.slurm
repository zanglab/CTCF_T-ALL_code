#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP4V2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4V2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4X1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4Z1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP4Z2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP51A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP7A1 
