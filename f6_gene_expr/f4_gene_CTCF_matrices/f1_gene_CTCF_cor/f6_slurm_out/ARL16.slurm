#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARL16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL17A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL17B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL2-SNX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL2BP 
