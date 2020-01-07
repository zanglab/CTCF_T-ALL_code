#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP46.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP50 
