#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UGT1A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A8 
