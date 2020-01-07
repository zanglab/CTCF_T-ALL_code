#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM81A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM81A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM81B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83D 
