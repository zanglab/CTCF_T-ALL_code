#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX51.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX56 
