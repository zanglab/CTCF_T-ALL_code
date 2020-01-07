#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TFPI.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TFPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFPI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFRC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TG 
