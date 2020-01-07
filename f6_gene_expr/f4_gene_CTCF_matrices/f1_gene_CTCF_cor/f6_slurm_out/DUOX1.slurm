#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUOX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUOXA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUOXA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUS1L 
