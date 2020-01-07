#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LETM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LETM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LETMD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LEUTX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LFNG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS12 
