#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AP4B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AP4B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP4E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP4M1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP4S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP5B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AP5M1 
