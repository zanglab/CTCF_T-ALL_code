#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLAC9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAGL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAGL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAU 
