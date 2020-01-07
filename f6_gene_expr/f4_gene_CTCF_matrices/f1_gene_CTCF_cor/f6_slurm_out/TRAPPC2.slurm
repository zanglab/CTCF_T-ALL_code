#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRAPPC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRAPPC4 
