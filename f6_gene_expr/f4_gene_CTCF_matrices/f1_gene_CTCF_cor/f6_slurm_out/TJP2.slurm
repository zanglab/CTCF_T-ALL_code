#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TJP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TJP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TJP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TKFC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TKT 
