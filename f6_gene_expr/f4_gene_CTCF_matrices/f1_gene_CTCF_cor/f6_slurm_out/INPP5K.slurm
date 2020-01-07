#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INPP5K.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INS-IGF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INSIG1 
