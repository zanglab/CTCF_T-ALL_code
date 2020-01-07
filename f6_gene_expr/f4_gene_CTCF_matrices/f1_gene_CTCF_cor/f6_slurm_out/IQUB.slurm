#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IQUB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IQUB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRAK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRAK1BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRAK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRAK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRAK4 
