#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPN12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN20 
