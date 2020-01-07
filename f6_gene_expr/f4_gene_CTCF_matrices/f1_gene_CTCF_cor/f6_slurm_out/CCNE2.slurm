#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCNE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNI 
