#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LTF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LTF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LUC7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LUC7L2 
