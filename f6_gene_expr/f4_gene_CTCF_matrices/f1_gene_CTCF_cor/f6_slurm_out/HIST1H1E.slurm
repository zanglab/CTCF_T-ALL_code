#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H1E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H1T 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AD 
