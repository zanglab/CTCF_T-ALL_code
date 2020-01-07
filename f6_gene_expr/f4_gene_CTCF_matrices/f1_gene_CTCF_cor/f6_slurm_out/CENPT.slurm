#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CENPT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPW 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP104 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP112 
