#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRCAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRCIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRD5A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRD5A1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRD5A2 
