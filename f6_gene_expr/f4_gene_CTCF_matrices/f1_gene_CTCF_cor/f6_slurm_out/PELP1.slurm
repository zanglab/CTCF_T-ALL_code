#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PELP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PELP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PENK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEPD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PER2 
