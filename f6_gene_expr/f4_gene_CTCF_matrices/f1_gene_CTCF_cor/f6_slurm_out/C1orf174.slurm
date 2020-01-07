#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1orf174.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf174 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf185 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf189 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf194 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf198 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf21 
