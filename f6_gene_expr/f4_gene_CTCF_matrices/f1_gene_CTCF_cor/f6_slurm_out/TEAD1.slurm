#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEAD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEAD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEAD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TECPR1 
