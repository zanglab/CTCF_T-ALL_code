#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OTX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OTX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OVCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OVCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OVCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OVGP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OVOL1 
