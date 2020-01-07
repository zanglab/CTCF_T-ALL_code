#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF93.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF93 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF98 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF99 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNFX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNHIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNHIT2 
