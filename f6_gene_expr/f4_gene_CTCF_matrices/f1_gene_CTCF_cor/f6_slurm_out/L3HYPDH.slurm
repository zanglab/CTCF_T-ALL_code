#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o L3HYPDH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i L3HYPDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L3MBTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L3MBTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L3MBTL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L3MBTL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LACC1 
