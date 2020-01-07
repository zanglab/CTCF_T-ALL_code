#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMPRSS11F.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11GP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS2 
