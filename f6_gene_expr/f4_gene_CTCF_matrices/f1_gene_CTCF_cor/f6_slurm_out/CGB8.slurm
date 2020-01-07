#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CGB8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CGB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGGBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGREF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CGRRF1 
