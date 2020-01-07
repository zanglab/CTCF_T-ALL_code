#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2V1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2V1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2V2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2W 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2Z 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE3B 
