#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZCCHC9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCWPW1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCWPW2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDBF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZDHHC1 
