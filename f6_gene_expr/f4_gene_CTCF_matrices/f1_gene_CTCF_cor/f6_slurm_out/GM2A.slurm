#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GM2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GM2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMEB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMEB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMFB 
