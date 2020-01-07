#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCL7C 
