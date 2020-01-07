#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FUT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT9 
