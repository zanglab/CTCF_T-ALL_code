#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PGLYRP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PGLYRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGLYRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGLYRP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGLYRP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGM2 
