#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM216A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM216A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM216B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM217A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM217B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM218A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM219A 
