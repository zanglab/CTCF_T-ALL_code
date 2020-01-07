#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSGA10IP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSGA10IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSGA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSHZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSHZ2 
