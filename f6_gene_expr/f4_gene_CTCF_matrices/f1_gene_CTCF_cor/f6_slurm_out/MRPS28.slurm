#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPS28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS31P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS34 
