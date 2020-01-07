#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLCC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLCN4 
