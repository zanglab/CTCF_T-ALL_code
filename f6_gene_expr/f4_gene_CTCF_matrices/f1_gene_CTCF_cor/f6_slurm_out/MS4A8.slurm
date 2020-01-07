#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MS4A8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSANTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSANTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSANTD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSANTD3-TMEFF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSANTD4 
