#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLRD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRK1 
