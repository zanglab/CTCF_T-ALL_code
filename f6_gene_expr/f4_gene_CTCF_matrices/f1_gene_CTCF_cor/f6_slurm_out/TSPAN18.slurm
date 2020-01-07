#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSPAN18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSPAN32 
