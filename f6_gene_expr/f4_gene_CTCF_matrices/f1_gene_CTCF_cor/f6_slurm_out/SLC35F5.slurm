#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC35F5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35F5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35F6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35G1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35G5 
