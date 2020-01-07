#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o H2AFY2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2BFM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2BFWT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H3F3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H3F3AP4 
