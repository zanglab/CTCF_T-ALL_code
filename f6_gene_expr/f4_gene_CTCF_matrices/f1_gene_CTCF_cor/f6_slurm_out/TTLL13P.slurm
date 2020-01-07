#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTLL13P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL13P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL6 
