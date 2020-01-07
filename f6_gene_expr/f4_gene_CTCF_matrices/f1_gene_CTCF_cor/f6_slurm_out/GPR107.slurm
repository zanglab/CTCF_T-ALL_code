#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR107.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR107 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR108 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR119 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR132 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR135 
