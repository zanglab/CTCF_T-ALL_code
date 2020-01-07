#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGMO.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGMO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGPAT1 
