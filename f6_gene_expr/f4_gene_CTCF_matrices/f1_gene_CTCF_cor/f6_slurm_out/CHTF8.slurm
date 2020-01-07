#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHTF8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHTF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHTOP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHUK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHURC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHURC1-FNTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIAO1 
