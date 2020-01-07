#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPRN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPRXL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPSAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPSB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPSD1 
