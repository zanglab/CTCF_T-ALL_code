#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHDC7A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHDC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL1 
