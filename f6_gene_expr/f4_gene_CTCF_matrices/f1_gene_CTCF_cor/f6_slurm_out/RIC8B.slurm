#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RIC8B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RIC8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RICTOR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RIIAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RILP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RILPL1 
