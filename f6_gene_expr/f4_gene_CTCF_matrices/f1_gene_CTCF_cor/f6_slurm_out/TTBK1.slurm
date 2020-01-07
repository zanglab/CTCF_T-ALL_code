#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTBK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTBK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTBK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC14 
