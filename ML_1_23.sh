### ML analysis ###

#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -N ML_1_23
#PBS -M 21634076@sun.ac.za  
#PBS -m abe
#PBS -e ML_1_23.err
#PBS -o ML_1_23.out
#PBS -l walltime=100:00:00

#Declare working directory from which to submit PBS script
cd $PBS_O_WORKDIR

# Load IQTree 1.6.12
module load app/IQTREE/1.6.12

iqtree -s /home/21634076/Carcharhiniform_mitogenomics/ML/Galeomorphi_cleaned_alignment.phy -spp /home/21634076/Carcharhiniform_mitogenomics/ML/partitions_1_23/partitions_1_23.NEX.best_scheme.nex -alrt 1000 -bb 1000