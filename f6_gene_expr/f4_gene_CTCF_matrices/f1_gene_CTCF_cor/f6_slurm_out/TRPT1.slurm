#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRPT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPV2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPV3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPV4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPV5 
