#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUSP27.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP5P1 
