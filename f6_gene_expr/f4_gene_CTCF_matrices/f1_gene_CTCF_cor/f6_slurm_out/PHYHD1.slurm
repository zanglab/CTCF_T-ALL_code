#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHYHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHYHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHYHIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHYHIPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHYKPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PI16 
