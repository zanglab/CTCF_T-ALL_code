#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MLIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MLIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLKL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLLT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLLT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLLT10P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLLT11 
