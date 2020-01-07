#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FRMD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD8P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMPD1 
