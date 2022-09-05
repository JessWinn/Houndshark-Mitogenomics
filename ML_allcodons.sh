### ML analysis ###

#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -N ML_allcodons
#PBS -M 21634076@sun.ac.za  
#PBS -m abe
#PBS -e ML_allcodons.err
#PBS -o ML_allcodons.out
#PBS -l walltime=100:00:00

#Declare working directory from which to submit PBS script
cd $PBS_O_WORKDIR

# Load IQTree 1.6.12
module load app/IQTREE/1.6.12

iqtree -s /home/21634076/Carcharhiniform_mitogenomics/ML/Galeomorphi_cleaned_alignment.phy -spp /home/21634076/Carcharhiniform_mitogenomics/ML/partitions_allcodons/partitions_allcodons.nex.best_scheme.nex -alrt 1000 -bb 1000