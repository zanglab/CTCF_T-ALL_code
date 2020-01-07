#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COG5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL10A1 
