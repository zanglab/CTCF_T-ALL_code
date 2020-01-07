#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNA10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA6 
