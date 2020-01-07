#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIX6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKAP2 
