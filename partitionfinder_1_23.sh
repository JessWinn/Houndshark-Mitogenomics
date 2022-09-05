### PARTITIONING & MODEL SELECTION ###

#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -N partitions_1_23
#PBS -M 21634076@sun.ac.za  
#PBS -m abe
#PBS -e partitions_1_23.err
#PBS -o partitions_1_23.out
#PBS -l walltime=01:00:00

#Declare working directory from which to submit PBS script
cd $PBS_O_WORKDIR

# Load IQTree 1.6.12
module load app/IQTREE/1.6.12

## 1. Export the final concatenated cleaned alignment as a phylip file from Geneious.

## 2. Create a NEXUS file containing annotation information

iqtree -s /home/21634076/Carcharhiniform_mitogenomics/ML/Galeomorphi_cleaned_alignment.phy -spp /home/21634076/Carcharhiniform_mitogenomics/ML/partitions_1_23/partitions_1_23.NEX -m MF+MERGE -rcluster 30 -AICc -safe -nt AUTO

