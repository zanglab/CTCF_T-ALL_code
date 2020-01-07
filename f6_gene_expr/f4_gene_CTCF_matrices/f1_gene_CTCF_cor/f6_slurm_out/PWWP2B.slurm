#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PWWP2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PWWP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXDNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXMP2 
