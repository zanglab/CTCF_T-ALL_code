#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TIPARP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TIPARP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIPIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIPRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIRAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TJAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TJP1 
