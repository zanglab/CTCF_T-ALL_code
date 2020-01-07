#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOP3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOPAZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOPBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOPORS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOR1AIP1 
