#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o B3GNT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT9 
