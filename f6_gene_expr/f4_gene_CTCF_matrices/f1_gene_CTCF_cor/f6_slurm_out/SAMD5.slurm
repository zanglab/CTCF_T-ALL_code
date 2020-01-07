#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SAMD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMD9L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SAMHD1 
