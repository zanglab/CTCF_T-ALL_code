#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DBH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DBH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBIL5P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBIL5P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBNDD1 
