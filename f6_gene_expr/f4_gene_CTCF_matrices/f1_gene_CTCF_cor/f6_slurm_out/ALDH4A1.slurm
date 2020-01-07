#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALDH4A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH4A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH5A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH6A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH7A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH8A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH9A1 
