#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUTM2F.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTM2F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUTM2G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NVL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NWD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NWD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXF1 
