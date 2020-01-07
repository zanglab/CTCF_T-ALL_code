#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CMTR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMTR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CMYA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNBP 
