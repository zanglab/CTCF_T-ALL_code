#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LDLRAD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LDLRAD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDLRAD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDLRAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LDOC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LECT2 
