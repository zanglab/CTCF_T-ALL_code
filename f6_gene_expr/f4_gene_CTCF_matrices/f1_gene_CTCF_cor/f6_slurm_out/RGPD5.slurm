#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGPD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS10 
