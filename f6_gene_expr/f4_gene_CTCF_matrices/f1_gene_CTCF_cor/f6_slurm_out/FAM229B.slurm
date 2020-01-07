#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM229B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM229B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM231C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM24A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM24B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM25A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM25BP 
