#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COL20A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COL20A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL21A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL22A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL23A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL24A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL25A1 
