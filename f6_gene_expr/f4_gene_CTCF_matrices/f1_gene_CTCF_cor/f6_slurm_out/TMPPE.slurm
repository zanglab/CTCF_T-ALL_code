#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMPPE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11BNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS11E 
