#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIR3DX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIR3DX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIRREL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIRREL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KISS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KISS1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIT 
