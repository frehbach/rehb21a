#!/bin/bash
#SBATCH --job-name=21a
#SBATCH --time=2-12:00:00 #2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=all #testing
#SBATCH --output=slurmOut/%A-%a.out
#SBATCH --exclude=s7cfd16
#SBATCH --mail-user=frederik.re@gmail.com
#SBATCH --mail-type=NONE

echo ~
source ~/.bashrc
echo $PATH
hostname
which R
