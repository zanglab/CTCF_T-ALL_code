#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP27X 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP3 
