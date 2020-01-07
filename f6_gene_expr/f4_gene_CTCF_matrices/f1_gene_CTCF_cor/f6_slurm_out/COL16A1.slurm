#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COL16A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COL16A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL17A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL18A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL19A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL1A2 
