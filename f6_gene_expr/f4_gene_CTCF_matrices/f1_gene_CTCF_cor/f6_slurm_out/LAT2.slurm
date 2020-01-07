#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LAT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LATS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LATS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAYN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LBH 
