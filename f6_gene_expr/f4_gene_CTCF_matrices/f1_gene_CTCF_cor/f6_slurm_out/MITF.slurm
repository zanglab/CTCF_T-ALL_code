#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MITF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MITF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIXL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MKI67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MKKS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MKL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MKL2 
