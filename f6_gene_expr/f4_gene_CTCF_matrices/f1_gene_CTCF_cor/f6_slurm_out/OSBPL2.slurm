#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OSBPL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL8 
