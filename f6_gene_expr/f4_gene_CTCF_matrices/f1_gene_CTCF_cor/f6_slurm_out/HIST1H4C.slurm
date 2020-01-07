#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H4C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4H 
