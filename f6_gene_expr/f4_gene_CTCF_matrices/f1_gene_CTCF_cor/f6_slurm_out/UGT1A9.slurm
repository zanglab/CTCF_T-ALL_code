#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UGT1A9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B11 
