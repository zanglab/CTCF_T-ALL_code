#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PEX13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX26 
