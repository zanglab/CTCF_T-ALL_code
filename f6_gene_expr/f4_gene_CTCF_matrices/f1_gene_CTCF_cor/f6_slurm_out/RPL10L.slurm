#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL10L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL10L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13AP17 
