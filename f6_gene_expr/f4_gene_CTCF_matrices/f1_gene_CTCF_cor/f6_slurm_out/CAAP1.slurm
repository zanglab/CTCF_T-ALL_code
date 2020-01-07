#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAB39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAB39L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABLES1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABLES2 
