#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM182B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM182B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM183A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM184A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM184B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM185A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM186A 
