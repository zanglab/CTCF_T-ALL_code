#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COPG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COPG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPRS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPS5 
