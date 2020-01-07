#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZBED6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBED6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBED6CL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBED8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBED9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB1 
