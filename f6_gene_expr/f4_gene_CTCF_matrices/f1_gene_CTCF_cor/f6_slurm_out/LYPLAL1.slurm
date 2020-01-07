#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LYPLAL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPLAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYRM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYRM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYRM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYRM9 
