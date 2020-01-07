#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYLK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYLK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYLPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYNN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO15A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO15B 
