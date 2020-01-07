#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AFMID.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AFMID 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AFTPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGAP2 
