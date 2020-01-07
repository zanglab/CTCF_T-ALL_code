#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCBP2-AS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCCRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCEH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCF1 
