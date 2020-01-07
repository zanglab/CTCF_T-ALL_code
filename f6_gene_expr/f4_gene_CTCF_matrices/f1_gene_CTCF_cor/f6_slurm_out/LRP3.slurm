#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP8 
