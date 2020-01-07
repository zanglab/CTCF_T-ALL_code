#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o F8A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i F8A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F8A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FA2H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAAH2 
