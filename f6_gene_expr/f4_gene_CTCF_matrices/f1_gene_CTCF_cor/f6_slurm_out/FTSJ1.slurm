#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FTSJ1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FTSJ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FTSJ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUCA2 
