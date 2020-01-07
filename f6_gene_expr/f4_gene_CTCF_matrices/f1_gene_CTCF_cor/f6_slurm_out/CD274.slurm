#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD274.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD274 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD276 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD2AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD2BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300A 
