#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FYCO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FYCO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FYN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FYTTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FZD2 
