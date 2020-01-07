#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM135B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM135B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM136A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM13C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM149A 
