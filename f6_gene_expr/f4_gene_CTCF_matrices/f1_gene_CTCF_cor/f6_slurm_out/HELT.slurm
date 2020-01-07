#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HELT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HELT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HELZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HELZ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEMGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEMK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HENMT1 
