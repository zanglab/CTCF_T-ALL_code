#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTER.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTER 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGDR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGDR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGER1 
