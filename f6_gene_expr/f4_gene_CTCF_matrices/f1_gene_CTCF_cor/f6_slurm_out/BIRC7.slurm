#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BIRC7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIVM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIVM-ERCC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLCAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BLID 
