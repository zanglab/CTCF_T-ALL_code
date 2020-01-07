#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab 
#SBATCH -o run2.out

time python py2_methy_compr.py
