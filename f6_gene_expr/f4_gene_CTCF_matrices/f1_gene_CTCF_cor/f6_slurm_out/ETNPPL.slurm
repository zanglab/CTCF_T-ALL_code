#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ETNPPL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ETNPPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV3 
