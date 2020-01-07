#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPOCD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPOCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPOCK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPOCK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPOCK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPON1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPON2 
