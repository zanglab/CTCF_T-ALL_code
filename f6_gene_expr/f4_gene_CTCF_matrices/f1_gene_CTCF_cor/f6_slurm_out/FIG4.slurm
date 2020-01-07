#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FIG4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FIG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIGLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIGNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIGNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FILIP1 
