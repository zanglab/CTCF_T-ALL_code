#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERRFI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERRFI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERV3-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERVFRD-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERVMER34-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERVV-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERVV-2 
