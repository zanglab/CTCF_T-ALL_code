#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTPS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTR9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTRC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTRL 
