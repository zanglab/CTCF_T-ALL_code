#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RRN3P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RRN3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRN3P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRNAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RRP15 
