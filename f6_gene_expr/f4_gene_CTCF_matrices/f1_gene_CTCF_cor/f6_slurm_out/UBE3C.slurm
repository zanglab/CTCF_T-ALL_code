#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE3C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBFD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBIAD1 
