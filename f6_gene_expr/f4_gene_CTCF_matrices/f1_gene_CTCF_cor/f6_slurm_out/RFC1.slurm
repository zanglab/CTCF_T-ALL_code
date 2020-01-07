#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RFC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RFC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RFESD 
