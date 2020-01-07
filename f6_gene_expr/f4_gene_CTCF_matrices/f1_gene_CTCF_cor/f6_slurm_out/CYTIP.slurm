#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYTIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYTIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYYR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i D2HGDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAAM2 
