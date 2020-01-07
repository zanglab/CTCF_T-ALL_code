#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSPYL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPYL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPYL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPYL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSR3 
