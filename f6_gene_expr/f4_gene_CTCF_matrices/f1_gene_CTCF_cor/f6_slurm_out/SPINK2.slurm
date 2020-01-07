#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPINK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK8 
