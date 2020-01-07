#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ICK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ICK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICOS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ICOSLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ID1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ID2 
