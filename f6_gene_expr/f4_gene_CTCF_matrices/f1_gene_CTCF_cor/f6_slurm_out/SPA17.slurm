#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPA17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPA17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPACA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPACA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPACA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPACA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPACA5B 
