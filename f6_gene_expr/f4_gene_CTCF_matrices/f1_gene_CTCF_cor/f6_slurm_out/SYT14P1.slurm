#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYT14P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT14P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT3 
