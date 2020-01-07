#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP34.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP39 
