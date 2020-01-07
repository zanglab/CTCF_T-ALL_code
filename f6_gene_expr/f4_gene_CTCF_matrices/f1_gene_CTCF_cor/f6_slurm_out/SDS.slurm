#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SDS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEBOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC11C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC13 
