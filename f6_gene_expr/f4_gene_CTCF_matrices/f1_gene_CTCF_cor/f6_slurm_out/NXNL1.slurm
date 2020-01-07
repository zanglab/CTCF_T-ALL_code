#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NXNL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NXNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPE4 
