#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NLGN4Y.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NLGN4Y 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRC5 
