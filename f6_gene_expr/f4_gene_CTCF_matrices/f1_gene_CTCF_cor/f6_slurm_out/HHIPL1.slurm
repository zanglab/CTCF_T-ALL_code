#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HHIPL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HHIPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHIPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHLA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHLA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HHLA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIBADH 
