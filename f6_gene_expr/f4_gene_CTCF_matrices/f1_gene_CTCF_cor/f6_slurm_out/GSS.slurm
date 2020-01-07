#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GSS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GSS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTA5 
