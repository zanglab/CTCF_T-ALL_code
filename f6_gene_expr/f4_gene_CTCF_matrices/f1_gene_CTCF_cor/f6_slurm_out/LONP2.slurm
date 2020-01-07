#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LONP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LONP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LONRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LONRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LONRF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LOR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LOX 
