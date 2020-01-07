#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTF3C5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3C5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF3C6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTPBP3 
