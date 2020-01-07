#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLK3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLLU1OS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLMN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLN3 
