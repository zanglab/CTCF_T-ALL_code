#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DZANK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DZANK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DZIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DZIP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DZIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i E2F2 
