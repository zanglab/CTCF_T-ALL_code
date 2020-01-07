#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DLGAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DLGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLGAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLGAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLGAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLGAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLK1 
