#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C16orf96.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf96 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf107 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C17orf50 
