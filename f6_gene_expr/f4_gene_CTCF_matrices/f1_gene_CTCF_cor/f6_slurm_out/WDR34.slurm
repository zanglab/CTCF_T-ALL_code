#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR34.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR4 
