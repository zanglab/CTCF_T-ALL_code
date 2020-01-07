#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COX16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COX16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX4I1 
