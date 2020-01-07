#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DES.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DESI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DESI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DET1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEXI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DFFA 
