#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX39A 
