#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR78.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR83 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR83OS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR86 
