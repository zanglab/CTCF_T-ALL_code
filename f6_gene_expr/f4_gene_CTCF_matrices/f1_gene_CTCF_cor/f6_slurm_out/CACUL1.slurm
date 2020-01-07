#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CACUL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CACUL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACYBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CADM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CADM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CADM3 
