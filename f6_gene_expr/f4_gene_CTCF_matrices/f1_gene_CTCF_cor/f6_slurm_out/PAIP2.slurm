#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAIP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAIP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAK1IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAK3 
