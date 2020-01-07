#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERCC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC6L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC6L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EREG 
