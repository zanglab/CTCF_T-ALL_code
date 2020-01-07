#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BTNL9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BTNL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTRC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BUB1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BUB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BUD13 
