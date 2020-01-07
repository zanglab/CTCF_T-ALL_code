#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAGE12F.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE12F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE12J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAGE2E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAK 
