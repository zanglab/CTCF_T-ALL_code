#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAPGEF4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPGEFL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAPSN 
