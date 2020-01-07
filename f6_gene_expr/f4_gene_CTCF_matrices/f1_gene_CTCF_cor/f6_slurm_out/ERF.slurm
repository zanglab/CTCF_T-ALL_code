#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERGIC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERGIC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERGIC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERH 
