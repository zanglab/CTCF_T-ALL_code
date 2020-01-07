#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HECW1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HECW1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HECW2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HELB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HELLS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HELQ 
