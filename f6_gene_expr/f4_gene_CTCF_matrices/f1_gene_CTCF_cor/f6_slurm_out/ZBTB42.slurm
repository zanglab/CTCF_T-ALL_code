#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZBTB42.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZBTB47 
