#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HAUS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS7 
