#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPN21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN5 
