#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PKIA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PKIA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKIB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKIG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKLR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PKMYT1 
