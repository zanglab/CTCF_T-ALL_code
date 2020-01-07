#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STOML3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STOML3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STON1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STON1-GTF2A1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STON2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STOX2 
