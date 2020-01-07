#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DGKE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DGKQ 
