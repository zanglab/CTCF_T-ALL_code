#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MGRN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MGRN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGST3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIA-RAB4B 
