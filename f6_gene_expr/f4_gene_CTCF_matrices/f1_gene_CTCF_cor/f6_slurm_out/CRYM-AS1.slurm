#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRYM-AS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYM-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYZL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSAG1 
