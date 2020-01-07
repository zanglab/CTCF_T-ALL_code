#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PER3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PERM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PERP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PES1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PET100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PET117 
