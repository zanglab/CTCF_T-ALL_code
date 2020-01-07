#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD300C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300LB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300LD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300LF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD300LG 
