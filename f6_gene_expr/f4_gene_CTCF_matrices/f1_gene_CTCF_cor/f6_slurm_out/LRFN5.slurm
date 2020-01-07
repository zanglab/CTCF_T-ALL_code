#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRFN5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRFN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRGUK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRIG2 
