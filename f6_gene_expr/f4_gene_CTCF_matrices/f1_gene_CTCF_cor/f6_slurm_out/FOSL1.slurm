#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FOSL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FOSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOSL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FOXB1 
