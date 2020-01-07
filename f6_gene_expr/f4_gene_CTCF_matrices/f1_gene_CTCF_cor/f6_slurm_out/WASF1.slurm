#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WASF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WASF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WASF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WASF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WASH2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WASH3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WASH5P 
