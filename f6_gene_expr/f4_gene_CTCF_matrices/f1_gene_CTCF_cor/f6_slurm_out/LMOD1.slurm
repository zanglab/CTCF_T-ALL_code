#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LMOD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LMOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMOD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMOD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMTK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMTK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMX1A 
