#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELAC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELANE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELAVL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELAVL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELAVL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELAVL4 
