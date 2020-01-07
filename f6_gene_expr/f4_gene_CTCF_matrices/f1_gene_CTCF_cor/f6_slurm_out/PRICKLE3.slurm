#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRICKLE3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRICKLE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRICKLE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRIMA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRIMPOL 
