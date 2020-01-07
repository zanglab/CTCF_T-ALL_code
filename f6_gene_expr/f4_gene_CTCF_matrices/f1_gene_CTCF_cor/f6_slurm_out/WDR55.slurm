#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR55.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR61 
