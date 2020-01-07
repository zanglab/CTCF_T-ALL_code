#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHL10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL15 
