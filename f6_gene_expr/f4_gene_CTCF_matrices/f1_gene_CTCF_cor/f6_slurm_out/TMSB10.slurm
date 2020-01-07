#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMSB10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMSB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMSB15A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMSB15B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMSB4X 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMSB4Y 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMTC1 
