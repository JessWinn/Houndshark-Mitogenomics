#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -N spades_de_novo
#PBS -M 21634076@sun.ac.za  
#PBS -m abe
#PBS -e spades_de_novo.err
#PBS -o spades_de_novo.out
#PBS -l walltime=900:00:00

#Declare working directory from which to submit PBS script
cd $PBS_O_WORKDIR

work=/home/21634076/Iontorrent/Spades_de_novo

#Load spades and quast
module load app/SPAdes
module load app/QUAST

#Mustelus palumbes
spades.py -k 21,33,55,77,99,127 -t 8 -m 250  --iontorrent -s $work/IonCode_0192_rawlib.basecaller.bam --careful -o $work/Mpa 
quast.py -o $work/Mpa/QuastMpa $work/Mpa/contigs.fasta 

#Mustelus mosis
spades.py -k 21,33,55,77,99,127 -t 8 -m 250  --iontorrent -s $work/IonCode_0193_rawlib.basecaller.bam --careful -o $work/Mmos 
quast.py -o $work/Mmos/QuastMmos $work/Mmos/contigs.fasta

#Triakis megalopterus
spades.py -k 21,33,55,77,99,127 -t 8 -m 250  --iontorrent -s $work/IonCode_0191_rawlib.basecaller.bam --careful -o $work/Tmeg 
quast.py -o $work/Tmeg/QuastTmeg $work/Tmeg/contigs.fasta

#Galeorhinus galeus
spades.py -k 21,33,55,77,99,127 -t 8 -m 250  --iontorrent -s $work/IonCode_0450_rawlib.basecaller.bam --careful -o $work/GG 
quast.py -o $work/GG/QuastGG $work/GG/contigs.fasta

#Mustelus asterias
spades.py -k 21,33,55,77,99,127 -t 8 -m 250  --iontorrent -s $work/IonCode_0449_rawlib.basecaller.bam --careful -o $work/Mast 
quast.py -o $work/Mast/QuastMast $work/Mast/contigs.fasta
