#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PYGO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PYGO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYGO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYHIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYROXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYROXD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYURF 
