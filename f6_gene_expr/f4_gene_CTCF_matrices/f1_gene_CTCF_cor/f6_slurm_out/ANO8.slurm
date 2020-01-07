#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANO8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANO9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANOS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANP32A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANP32AP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANP32B 
