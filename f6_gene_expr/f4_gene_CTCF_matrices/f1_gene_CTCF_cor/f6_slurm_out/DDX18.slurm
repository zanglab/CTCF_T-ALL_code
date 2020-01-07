#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX19A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX19B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX23 
