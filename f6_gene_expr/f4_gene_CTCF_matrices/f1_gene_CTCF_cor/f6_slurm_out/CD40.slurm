#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD40LG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD48 
