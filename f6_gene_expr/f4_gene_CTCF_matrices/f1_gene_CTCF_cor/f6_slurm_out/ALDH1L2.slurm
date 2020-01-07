#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALDH1L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH3A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH3A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH3B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH3B2 
