#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MNAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MNAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MNDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MNS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MNT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MNX1 
