#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPEB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPEB2 
