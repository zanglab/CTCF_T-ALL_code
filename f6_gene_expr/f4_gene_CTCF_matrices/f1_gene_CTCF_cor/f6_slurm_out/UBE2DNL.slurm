#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2DNL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2DNL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2E3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2F-SCLY 
