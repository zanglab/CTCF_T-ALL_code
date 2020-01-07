#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UGT2B15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT2B7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT3A1 
