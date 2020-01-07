#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP9X 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP9Y 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTF1 
