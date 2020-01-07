#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLTPD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLTPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLUD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLUD1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLUD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLUL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYAT 
