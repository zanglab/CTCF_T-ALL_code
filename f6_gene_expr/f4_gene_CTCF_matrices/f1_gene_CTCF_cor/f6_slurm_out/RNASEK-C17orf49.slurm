#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNASEK-C17orf49.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEK-C17orf49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASET2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RND2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RND3 
