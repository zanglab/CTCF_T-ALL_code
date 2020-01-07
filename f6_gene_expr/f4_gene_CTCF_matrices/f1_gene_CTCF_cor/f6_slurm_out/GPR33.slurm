#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR37L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR39 
