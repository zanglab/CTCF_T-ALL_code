#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o P3H4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i P3H4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P4HA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P4HA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P4HA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P4HB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P4HTM 
