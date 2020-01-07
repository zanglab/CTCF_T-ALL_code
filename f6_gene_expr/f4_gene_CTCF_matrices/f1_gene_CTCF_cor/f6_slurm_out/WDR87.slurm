#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR87.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR89 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR90 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR92 
