#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SHISA8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHISA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHKBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHOC2 
