#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATF5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF7IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATF7IP2 
