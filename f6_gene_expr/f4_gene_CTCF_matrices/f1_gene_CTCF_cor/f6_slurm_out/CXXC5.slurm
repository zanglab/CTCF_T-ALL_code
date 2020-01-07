#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CXXC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CXXC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf40A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXorf40B 
