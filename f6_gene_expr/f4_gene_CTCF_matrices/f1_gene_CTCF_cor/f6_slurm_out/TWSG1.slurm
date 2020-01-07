#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TWSG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TWSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXLNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXLNB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXLNG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TXLNGY 
