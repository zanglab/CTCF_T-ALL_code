#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR93.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR93 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDSUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDTC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDYHV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WEE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WEE2 
