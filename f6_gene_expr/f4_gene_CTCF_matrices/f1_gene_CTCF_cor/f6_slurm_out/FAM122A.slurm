#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM122A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM122A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM122B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM122C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM124A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM124B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM126A 
