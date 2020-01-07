#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSPY26P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPY26P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPY8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPYL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPYL2 
