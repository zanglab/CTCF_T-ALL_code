#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFAND4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAND6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFC3H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFHX2 
