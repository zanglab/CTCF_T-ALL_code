#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOLGA8N.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8N 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8O 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8S 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGA8T 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOLGB1 
