#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR33 
