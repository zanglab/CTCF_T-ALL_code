#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST2H2AB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2AB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2AC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2BA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2BC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2BE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2BF 
