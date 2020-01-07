#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP10B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP10D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP11C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP12A 
