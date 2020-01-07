#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP51.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP6NL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP7 
