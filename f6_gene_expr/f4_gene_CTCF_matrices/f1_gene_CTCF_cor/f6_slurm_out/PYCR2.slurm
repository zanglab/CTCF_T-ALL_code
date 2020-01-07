#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PYCR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PYCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYGM 
