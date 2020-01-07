#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C22orf23.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C22orf23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C22orf31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C22orf39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C22orf42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C22orf46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2CD2 
