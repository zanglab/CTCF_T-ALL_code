#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COX4I2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COX4I2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX6A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX6A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COX6B1 
