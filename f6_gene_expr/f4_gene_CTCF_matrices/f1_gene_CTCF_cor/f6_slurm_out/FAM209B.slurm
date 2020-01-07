#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM209B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM209B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM20A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM20B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM20C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM210A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM210B 
