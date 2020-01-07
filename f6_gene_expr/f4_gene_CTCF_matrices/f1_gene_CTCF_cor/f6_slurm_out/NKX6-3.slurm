#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NKX6-3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NKX6-3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLGN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLGN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLGN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLGN4X 
