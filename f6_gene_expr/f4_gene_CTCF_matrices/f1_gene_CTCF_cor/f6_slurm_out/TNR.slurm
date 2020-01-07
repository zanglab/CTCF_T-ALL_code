#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNRC18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNRC18P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNRC6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNRC6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNRC6C 
