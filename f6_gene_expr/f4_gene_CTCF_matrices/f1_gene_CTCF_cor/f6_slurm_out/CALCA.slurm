#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CALCA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CALCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALCOCO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALCOCO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALCR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALCRL 
