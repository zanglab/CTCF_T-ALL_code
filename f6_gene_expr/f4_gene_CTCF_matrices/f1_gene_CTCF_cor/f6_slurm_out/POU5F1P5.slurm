#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POU5F1P5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POU5F1P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU5F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU6F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POU6F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PP2D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPA1 
