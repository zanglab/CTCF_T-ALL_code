#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR151.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR151 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR152 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR153 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR155 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR156 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR157 
