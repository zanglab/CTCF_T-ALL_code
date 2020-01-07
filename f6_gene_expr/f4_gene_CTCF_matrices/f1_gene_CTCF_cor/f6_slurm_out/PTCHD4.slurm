#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTCHD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCHD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTDSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTDSS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTEN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTENP1 
