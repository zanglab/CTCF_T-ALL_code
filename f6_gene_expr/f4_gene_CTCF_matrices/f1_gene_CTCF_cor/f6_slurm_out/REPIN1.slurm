#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o REPIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i REPIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RERE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RERG 
