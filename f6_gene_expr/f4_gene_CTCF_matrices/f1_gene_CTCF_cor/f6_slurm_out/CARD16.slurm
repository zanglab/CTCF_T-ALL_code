#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CARD16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD9 
