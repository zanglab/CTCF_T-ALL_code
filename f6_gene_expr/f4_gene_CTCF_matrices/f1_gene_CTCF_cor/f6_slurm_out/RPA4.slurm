#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPA4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPAIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPE 
