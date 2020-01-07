#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1orf115.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf115 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf116 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf122 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf123 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf127 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf131 
