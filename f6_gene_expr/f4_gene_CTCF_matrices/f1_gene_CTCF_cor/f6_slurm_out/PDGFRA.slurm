#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDGFRA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDGFRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDHA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDHA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDHB 
