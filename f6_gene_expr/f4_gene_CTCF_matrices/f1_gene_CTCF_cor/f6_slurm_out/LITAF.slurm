#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LITAF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LITAF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIX1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LKAAEAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LLGL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LLGL2 
