#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LCMT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LCMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCMT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN15 
