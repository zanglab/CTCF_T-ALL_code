#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM92A1P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM92A1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM92B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM95C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM96A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM96B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM98A 
