#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCUN1D5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCUN1D5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCXR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDAH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDAH2 
