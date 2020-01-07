#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRA3 
