#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NXPH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXPH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NXT2 
