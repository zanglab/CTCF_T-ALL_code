#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTPBP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTSE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTSF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTSF1L 
