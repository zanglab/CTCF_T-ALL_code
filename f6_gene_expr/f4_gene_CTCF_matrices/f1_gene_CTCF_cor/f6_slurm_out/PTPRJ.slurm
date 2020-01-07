#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPRJ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRO 
