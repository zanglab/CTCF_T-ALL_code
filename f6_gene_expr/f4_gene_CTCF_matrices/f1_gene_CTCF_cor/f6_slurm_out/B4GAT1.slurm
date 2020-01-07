#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o B4GAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i B4GAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B9D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B9D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAALC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BABAM1 
