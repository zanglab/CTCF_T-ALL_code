#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAHD2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAHD2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAHD2CP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAIM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM102A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM102B 
