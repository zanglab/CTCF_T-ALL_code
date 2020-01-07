#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYOT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOZ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYOZ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYPOP 
