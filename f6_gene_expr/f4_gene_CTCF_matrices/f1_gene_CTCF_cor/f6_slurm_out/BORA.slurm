#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BORA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BORA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BORCS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BORCS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BORCS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BORCS7-ASMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BORCS8 
