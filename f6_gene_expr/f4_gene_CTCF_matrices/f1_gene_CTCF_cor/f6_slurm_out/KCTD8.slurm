#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCTD8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDELC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDELC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDELR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDELR2 
