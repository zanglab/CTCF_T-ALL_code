#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNPEP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNPEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNPEPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROBO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROBO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROBO3 
