#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCAMP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAMP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAND2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCAPER 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCARA3 
