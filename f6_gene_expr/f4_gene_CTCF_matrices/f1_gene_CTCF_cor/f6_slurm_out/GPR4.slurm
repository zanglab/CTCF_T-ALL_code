#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR6 
