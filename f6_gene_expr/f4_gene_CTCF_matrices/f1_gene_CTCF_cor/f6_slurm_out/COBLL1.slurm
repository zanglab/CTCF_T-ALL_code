#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COBLL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COBLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COG4 
