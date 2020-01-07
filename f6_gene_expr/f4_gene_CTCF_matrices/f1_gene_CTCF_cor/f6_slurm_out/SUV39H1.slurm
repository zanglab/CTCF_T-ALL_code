#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUV39H1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUV39H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUV39H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUZ12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUZ12P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SV2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SV2B 
