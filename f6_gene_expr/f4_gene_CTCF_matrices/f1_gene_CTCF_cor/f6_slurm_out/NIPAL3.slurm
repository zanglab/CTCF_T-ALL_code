#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NIPAL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPAL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPSNAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPSNAP3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NIPSNAP3B 
