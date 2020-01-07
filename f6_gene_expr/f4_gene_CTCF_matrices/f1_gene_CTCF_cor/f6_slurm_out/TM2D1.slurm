#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TM2D1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TM2D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM2D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM2D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TM4SF19 
