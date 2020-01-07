#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C2orf49.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf66 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf69 
