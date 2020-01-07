#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST2H3C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST3H2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST3H2BB 
