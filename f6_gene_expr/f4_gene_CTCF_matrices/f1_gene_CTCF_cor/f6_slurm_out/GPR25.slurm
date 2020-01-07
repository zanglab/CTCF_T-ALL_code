#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR32 
