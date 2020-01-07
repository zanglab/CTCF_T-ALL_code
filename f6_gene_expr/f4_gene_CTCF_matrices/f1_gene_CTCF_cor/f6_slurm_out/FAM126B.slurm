#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM126B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM126B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM129A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM129B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM129C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM131A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM131B 
