#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLP1R.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLP1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLP2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRA4 
