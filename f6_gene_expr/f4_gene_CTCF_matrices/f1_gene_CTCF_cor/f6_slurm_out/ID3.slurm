#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ID3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ID3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ID4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDH3A 
