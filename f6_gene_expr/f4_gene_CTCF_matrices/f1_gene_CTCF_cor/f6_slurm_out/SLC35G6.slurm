#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC35G6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35G6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC36A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC36A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC36A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC36A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC37A1 
