#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C10orf95.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf95 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf99 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C11orf24 
