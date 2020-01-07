#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COL11A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COL11A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL11A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL12A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL13A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL14A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COL15A1 
