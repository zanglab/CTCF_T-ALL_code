#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM155A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM155A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM155B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM156A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM156B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM157A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM157B 
