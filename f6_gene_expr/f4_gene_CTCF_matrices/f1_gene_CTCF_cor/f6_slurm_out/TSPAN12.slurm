#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSPAN12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN17 
