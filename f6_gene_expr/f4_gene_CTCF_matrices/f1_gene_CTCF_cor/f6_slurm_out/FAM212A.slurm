#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM212A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM212A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM212B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM213A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM213B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM214A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM214B 
