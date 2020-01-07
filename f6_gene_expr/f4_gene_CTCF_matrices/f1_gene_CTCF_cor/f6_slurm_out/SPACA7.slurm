#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPACA7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPACA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAG17 
