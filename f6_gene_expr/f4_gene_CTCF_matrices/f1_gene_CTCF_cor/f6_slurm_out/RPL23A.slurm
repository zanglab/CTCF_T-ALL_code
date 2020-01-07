#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL23A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23AP32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23AP53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23AP64 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23AP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL23AP82 
