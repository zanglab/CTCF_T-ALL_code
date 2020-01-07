#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DRD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DRD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRD5P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRGX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRICH1 
