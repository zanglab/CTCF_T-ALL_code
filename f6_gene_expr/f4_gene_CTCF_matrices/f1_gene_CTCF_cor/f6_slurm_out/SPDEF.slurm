#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPDEF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDEF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPDYE2 
