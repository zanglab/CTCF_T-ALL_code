#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAP2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAP2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEF3 
