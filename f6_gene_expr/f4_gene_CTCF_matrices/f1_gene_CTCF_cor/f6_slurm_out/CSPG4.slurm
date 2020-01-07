#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSPG4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSPG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSPG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSRNP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSRNP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSRNP3 
