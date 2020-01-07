#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM205C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM205C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM206A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM207A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM208A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM208B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM209A 
