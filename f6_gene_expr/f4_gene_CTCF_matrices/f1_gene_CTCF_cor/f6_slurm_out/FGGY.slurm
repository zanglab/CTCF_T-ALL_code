#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGGY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGGY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHAD1 
