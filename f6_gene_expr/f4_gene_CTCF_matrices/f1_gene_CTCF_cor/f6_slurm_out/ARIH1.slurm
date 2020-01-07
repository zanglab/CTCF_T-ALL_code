#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARIH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARIH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARIH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARIH2OS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL11 
