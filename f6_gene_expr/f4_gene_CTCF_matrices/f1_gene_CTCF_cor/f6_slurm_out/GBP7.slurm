#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GBP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GBP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCAT 
