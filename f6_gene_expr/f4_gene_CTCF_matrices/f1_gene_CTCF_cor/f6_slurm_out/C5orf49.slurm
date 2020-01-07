#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C5orf49.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf63 
