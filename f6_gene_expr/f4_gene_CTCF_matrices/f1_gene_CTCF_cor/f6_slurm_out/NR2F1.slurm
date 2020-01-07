#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NR2F1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2F6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR3C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR3C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR4A1 
