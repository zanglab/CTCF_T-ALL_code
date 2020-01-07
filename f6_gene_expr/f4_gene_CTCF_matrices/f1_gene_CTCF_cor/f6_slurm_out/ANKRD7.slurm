#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKS1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKS1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKS4B 
