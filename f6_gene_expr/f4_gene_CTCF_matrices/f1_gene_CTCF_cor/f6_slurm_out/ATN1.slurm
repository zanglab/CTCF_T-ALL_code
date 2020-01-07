#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATOH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATOH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATOH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP10A 
