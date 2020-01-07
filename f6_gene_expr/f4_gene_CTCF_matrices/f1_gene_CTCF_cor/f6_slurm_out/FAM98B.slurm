#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM98B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM98B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM98C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM9C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAN1 
