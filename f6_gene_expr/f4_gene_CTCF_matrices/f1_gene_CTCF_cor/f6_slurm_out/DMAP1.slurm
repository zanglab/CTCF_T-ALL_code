#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DMAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DMAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMBT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMBT1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMBX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMD 
