#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBQLN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBQLN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBQLN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBQLNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBR3 
