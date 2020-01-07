#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCLRE1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCLRE1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCLRE1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCLRE1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCP1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCP1B 
