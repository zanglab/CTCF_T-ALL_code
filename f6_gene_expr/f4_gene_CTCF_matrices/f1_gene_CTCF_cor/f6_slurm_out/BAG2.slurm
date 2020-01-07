#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BAG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BAG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAGE2 
