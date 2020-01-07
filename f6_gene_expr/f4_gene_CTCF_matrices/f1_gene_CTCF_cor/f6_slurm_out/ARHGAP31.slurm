#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGAP31.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP39 
