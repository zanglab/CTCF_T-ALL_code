#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NISCH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NISCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NKAIN3 
