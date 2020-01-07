#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGFG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGFG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGFG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGGF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGMAT 
