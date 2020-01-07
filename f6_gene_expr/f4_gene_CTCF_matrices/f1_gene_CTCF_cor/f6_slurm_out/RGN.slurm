#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGPD4 
