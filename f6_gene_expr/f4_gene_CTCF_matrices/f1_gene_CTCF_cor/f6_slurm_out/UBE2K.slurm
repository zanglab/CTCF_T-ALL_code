#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2K.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2MP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2N 
