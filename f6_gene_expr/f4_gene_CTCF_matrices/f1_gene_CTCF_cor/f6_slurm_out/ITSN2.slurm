#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITSN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITSN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IVD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IVL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IVNS1ABP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IWS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IYD 
