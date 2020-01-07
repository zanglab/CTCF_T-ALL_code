#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THUMPD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THUMPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THYN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TIAL1 
