#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAPT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAPVD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GARNL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GART 
