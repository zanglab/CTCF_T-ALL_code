#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o QRICH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i QRICH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QRSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QSER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QSOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QSOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QTRT1 
