#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USH1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USH1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USH2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USHBP1 
