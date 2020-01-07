#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o JMJD1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i JMJD1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JMJD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JMJD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JMJD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JMJD7-PLA2G4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JMJD8 
