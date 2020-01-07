#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CBWD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CBWD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBWD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBWD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBWD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBWD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBX1 
