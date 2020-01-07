#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMOT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMOTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMOTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMPD1 
