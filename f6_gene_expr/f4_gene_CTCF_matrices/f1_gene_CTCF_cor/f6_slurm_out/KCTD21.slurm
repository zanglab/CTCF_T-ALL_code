#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCTD21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD7 
