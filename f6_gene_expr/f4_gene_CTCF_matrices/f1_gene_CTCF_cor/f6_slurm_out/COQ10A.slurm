#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COQ10A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COQ10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COQ10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COQ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COQ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COQ4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COQ5 
