#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RRP1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP7BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP9 
