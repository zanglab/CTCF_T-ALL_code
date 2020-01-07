#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABLIM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABLIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABRACL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABT1 
