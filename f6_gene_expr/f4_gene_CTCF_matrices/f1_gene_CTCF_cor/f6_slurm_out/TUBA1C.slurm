#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TUBA1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBA4B 
