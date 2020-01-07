#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBC1D2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D22A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D22B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D23 
