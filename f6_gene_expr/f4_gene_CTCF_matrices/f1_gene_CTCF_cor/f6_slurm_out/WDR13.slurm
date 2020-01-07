#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR24 
