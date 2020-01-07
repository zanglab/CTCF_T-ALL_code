#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHAF1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD2 
