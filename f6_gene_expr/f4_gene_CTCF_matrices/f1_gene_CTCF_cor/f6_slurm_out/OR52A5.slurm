#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR52A5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52B6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR52E4 
