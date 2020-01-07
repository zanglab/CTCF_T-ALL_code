#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FSCB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FSCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSCN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FSD1L 
