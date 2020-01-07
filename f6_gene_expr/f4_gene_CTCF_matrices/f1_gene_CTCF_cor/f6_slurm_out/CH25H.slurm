#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CH25H.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CH25H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHADL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAF1A 
