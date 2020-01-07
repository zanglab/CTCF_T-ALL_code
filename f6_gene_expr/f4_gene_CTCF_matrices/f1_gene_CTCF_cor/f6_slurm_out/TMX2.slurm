#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMX2-CTNND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNF 
