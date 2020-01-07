#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CXorf49B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf49B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf67 
