#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C12orf49.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf60 
