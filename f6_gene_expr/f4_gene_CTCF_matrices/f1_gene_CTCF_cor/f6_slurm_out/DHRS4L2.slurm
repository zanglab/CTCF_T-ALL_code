#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DHRS4L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS4L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS7C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRS9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHRSX 
