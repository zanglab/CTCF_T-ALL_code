#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC35E3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35E3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35E4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35F1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35F2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35F3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC35F4 
