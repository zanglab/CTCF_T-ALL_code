#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHL35.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL41 
