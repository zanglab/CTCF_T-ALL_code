#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGFBP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBPL1 
