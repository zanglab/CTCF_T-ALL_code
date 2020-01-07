#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GMFG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GMFG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GML 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMNC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMNN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GMPPA 
