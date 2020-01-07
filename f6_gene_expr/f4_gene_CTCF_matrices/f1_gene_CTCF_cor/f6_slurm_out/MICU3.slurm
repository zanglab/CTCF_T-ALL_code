#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MICU3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MICU3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MID1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MID1IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MID2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIEF1 
