#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL17-C18orf32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL18A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL19 
