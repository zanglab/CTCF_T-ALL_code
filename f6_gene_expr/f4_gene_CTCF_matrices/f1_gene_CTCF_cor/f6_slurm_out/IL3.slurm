#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL31RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL34 
