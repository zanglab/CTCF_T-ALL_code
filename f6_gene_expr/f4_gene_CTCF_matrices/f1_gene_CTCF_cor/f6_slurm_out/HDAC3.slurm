#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HDAC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HDAC8 
