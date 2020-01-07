#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACTG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTG1P17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTG1P20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTG1P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTL10 
