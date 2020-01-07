#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C3orf22.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C3orf49 
