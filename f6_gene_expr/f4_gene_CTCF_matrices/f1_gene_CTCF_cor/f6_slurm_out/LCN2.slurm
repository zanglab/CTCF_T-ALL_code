#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LCN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCN9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCOR 
