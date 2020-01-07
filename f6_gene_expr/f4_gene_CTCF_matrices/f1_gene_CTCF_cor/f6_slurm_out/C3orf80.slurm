#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C3orf80.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf84 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4BPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4BPB 
