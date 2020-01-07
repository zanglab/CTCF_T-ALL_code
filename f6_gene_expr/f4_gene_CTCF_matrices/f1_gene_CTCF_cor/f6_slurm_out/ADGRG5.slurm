#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADGRG5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRL3 
