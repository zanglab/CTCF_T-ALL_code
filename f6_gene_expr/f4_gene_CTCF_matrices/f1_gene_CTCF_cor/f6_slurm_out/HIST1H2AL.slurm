#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H2AL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2APS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BC 
