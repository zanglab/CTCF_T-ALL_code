#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COPS6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPZ2 
