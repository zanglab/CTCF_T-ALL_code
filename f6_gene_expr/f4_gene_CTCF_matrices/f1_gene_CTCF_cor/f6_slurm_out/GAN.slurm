#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GANAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GANC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAP43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAPDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAPDHS 
