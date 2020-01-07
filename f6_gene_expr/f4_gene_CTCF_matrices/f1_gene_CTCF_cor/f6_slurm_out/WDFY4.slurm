#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDFY4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDFY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDPCP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR12 
