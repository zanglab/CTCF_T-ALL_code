#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDE9A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFD 
