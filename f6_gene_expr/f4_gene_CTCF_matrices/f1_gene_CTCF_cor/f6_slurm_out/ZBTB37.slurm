#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZBTB37.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB41 
