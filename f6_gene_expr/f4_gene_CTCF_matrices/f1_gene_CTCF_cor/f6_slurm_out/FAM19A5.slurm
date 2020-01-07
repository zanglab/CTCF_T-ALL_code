#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM19A5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM19A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM200A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM200B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM204A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM205A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM205BP 
