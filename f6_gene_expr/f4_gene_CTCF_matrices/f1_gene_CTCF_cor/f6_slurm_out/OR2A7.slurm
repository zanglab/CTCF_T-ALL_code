#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2A9P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2AE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2AG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2AG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2AK2 
