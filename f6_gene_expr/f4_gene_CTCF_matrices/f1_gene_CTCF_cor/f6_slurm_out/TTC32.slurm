#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTC32.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC38 
