#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPN9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRC 
