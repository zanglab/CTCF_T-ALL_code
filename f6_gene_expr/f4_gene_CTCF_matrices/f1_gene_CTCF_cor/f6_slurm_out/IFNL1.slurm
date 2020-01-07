#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFNL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNLR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNW1 
