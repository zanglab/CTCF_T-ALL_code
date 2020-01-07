#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BIN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BIRC6-AS2 
