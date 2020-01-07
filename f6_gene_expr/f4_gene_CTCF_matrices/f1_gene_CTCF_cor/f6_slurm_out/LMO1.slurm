#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LMO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LMO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMO7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMO7DN 
