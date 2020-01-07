#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL13AP20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13AP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13AP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13AP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13AP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL13P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL14 
