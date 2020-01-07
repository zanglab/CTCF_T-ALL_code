#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYB5D1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5R3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5R4 
