#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTP4A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTP4A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTP4A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN11 
