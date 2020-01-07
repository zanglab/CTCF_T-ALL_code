#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FMOD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FMOD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMR1NB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FN3K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FN3KRP 
