#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD70.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD79A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD79B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD80 
