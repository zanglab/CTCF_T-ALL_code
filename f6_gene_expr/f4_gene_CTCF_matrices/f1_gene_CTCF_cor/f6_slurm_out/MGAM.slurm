#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MGAM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGARP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT3 
