#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBXO42.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBXO47 
