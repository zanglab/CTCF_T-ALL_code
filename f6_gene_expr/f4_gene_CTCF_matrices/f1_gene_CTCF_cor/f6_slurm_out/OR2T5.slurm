#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR2T5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2T8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2V2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR2W3 
