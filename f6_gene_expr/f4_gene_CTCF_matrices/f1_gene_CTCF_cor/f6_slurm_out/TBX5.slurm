#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBX5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBXA2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBXAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TC2N 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCAF1 
