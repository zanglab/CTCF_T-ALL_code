#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLI.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLQ 
