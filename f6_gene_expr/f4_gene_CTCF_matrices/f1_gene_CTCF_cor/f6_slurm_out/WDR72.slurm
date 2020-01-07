#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WDR72.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR75 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WDR77 
