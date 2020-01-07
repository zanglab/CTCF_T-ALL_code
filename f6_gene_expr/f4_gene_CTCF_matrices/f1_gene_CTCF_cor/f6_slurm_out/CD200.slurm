#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD200.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD200 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD200R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD200R1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD207 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD209 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD22 
