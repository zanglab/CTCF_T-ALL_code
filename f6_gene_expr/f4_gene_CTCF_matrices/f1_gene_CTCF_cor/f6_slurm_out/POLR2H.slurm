#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLR2H.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2I 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2J2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2J3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2J4 
