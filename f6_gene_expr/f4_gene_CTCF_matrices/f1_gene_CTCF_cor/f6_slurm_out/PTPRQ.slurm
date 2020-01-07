#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPRQ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRVP 
