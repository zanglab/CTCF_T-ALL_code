#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRASP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRASP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRB14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRB7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GREB1 
