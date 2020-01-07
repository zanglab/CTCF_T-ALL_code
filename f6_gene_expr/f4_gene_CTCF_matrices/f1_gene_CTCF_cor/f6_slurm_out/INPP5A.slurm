#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INPP5A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INPP5J 
