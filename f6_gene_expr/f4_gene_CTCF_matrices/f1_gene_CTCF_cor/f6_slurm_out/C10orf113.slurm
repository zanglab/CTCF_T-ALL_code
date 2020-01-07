#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C10orf113.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf113 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf128 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf142 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf62 
