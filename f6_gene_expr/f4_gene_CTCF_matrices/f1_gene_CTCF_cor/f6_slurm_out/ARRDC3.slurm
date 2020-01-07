#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARRDC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRDC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARRDC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSD 
