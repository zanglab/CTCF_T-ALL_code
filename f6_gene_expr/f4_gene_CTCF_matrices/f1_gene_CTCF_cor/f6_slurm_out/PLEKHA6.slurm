#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLEKHA6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHA8P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLEKHB2 
