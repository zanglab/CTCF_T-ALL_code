#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR47.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR54 
