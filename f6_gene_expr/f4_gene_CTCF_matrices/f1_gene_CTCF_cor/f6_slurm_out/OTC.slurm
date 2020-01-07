#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OTC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OTC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOL1 
