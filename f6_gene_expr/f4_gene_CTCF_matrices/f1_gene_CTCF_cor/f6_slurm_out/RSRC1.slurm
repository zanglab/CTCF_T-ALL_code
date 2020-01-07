#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RSRC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RSRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSRC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSU1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTBDN 
