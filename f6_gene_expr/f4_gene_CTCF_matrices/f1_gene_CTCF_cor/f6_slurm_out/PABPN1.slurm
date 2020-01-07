#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PABPN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PABPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PABPN1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PACRG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PACRGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PACS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PACS2 
