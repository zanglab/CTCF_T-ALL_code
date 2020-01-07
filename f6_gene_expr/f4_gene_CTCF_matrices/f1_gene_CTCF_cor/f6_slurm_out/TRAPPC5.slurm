#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRAPPC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAT1 
