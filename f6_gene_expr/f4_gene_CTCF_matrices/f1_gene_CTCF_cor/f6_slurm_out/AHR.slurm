#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AHR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHRR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHSA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHSG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHSP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AICDA 
