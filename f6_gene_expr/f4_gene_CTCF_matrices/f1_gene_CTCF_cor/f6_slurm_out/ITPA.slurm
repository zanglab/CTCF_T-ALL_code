#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITPA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPKA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPKC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPR1 
