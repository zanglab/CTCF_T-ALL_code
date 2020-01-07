#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OSBPL9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OSBPL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSCAR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSGEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSGEPL1 
