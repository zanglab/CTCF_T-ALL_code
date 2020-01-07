#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR158.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR158 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR160 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR161 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR162 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR171 
