#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR62.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR64 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR66 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR70 
