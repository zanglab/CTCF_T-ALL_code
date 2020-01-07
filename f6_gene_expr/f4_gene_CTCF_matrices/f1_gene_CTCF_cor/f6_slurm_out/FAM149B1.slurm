#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM149B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM149B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM151A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM151B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM153A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM153B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM153C 
