#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKMY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKMY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKMY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD11 
