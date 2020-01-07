#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR7E47P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E47P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7E91P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR7G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR8B12 
