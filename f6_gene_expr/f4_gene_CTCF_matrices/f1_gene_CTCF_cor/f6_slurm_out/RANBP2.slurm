#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RANBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANGAP1 
