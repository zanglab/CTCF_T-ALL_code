#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NFKB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKBIA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKBIB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKBID 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKBIE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKBIL1 
