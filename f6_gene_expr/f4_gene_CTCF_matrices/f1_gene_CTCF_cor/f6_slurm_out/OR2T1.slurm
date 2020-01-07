#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2T1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T27 
