#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BAHCC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BAHCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAIAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAIAP2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAIAP2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAIAP3 
