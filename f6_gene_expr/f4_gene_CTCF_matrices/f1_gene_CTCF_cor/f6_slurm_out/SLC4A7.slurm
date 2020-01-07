#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC4A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC50A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC51A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC51B 
