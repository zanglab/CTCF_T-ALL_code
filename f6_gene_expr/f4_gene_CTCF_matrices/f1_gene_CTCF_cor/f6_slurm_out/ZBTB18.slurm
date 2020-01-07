#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZBTB18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB24 
