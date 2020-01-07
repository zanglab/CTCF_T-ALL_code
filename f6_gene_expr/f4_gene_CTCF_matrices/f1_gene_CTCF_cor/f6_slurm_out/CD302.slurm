#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD302.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD302 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD320 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD37 
