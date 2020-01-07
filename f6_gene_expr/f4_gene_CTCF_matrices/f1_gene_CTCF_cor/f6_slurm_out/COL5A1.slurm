#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COL5A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COL5A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL5A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL5A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A3 
