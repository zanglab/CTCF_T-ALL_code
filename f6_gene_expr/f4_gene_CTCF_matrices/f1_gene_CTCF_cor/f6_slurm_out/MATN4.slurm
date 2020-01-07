#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MATN4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MATN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MATR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAU2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAVS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAZ 
