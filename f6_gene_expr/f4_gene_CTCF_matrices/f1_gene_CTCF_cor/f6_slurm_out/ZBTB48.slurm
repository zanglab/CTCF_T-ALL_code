#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZBTB48.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB7B 
