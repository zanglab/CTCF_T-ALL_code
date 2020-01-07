#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPUSD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPUSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPUSD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPUSD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPUSD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRAGA 
