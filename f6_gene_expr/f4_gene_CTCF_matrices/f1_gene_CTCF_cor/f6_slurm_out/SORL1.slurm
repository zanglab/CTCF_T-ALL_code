#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SORL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SORL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SORT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOSTDC1 
