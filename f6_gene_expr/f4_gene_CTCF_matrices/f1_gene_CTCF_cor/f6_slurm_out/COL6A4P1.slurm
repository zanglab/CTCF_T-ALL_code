#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COL6A4P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A4P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A4P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL6A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL7A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL8A1 
