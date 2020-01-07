#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ODF3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF3L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ODF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OFD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGDH 
