#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCSK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCSK9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCTP 
