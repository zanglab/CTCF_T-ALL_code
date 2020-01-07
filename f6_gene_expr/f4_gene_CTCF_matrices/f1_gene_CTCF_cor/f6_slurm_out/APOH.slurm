#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APOH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APOH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APOL5 
