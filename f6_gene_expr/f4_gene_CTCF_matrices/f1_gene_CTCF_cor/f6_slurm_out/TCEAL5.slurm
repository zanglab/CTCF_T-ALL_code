#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCEAL5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEAL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEANC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCEANC2 
