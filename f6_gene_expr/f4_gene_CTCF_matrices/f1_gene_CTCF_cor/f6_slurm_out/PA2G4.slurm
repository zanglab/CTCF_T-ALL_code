#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PA2G4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PA2G4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PA2G4P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PABPC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PABPC1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PABPC1L2A 
