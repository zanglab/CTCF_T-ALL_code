#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR143.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR143 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR146 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR148 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR149 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR150 
