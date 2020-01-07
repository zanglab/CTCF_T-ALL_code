#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIBCH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIBCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HID1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIF1AN 
