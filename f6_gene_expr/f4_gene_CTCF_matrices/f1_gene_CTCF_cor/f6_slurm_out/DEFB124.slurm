#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DEFB124.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB124 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB125 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB126 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB132 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB133 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB134 
