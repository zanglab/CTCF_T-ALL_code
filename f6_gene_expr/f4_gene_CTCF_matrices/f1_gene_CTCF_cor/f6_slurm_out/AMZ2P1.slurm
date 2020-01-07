#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMZ2P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMZ2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC15 
