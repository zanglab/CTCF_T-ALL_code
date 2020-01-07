#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CASD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CASD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASKIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASKIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CASP10 
