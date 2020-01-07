#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR10H1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10H4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10H5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10J1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10J3 
