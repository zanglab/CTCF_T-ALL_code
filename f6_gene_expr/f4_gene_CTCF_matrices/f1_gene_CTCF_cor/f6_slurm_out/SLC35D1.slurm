#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC35D1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35E2B 
