#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FMO2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO6P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO9P 
