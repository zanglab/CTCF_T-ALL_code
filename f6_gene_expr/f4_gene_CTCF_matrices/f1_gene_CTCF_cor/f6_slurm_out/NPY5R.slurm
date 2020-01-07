#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPY5R.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPY5R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPY6R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NQO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NQO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR0B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR0B2 
