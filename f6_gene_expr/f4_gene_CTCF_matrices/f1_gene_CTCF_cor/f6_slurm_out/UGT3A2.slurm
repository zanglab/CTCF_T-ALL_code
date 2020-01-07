#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UGT3A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT3A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UHMK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UHRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UHRF1BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UHRF1BP1L 
