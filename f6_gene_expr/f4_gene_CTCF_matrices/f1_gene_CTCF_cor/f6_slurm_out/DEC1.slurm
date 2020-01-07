#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DEC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DEC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DECR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DECR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEDD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEDD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEF6 
