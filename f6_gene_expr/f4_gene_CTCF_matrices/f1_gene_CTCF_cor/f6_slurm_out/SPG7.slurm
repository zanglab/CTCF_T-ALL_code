#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPG7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPHAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPHK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPHK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPHKAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPI1 
