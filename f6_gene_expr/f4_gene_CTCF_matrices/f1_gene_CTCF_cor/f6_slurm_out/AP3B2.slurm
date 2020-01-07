#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AP3B2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AP3B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP3D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP3M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP3M2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP3S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP3S2 
