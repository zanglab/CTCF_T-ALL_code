#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPRCAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRCAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRH 
