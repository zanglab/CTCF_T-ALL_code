#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTCAP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTCAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTDAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KSR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KTI12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KTN1 
