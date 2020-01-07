#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COL26A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COL26A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL27A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL28A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL3A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL4A1 
