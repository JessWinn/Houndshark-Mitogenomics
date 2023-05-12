# Bioinformatics pipeline of Winn et al. (2022) - Triakid Mitophlyogenomics

## STEP 1: Quality Control of Ion GeneStudio™ S5 data

1. Check sequence quality in FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Trim adaptors and poor-quality bases (phred score below 16) and remove reads shorter than 25 base pairs (bp) with Torrent Suite Software 5.16.

## STEP 2: Mitogenome assemblies  

### STEP 2.1: reference-based assembly

1. Align raw reads to the _Mustelus mustelus_ mitogenome (NC_039629.1) using the Geneious read mapper with medium sensitivity settings and five iterations in Geneious Prime (version 2019.1.3).

### STEP 2.2: baited assembly

1. Extract reads that mapped to the reference genome in bam format.
2. Feed the reads into a de novo pipeline in SPAdes version 3.15 with the input set for unpaired Ion Torrent reads with 8 threads, kmers 21,33,55,77,99,127, the careful option to reduce the number of mismatches and short indels and all other parameters left as default.
EXAMPLE - Mustelus palumbes
```
#!/bin/bash
module load app/SPAdes
spades.py \
	-k 21,33,55,77,99,127 \
	-t 8 \
	-m 250 \  
	--iontorrent \ 
	-s $work/Mpa.bam \
	--careful \
	-o $work/Mpa \
```

### STEP 2.3: de novo assembly

Directly map raw reads (bam format) de novo.
EXAMPLE - Mustelus palumbes
```
#!/bin/bash
module load app/SPAdes
spades.py \
	-k 21,33,55,77,99,127 \
	-t 8 \
	-m 250 \  
	--iontorrent \ 
	-s $work/IonCode_0192_rawlib.basecaller.bam \
	--careful \
	-o $work/Mpa \
```

### STEP 2.4: Assembly comparison

1. Align the three assemblies to each other using the Geneious alignment tool with default parameters. 
2. Check each alignment for discrepancies on  Geneious and ediit maually to obtain the final genome sequence.

## STEP 3: Mitogenome annotations

1. Annotate protein coding genes (PCGs), ribosomal (r)RNA and transfer (t)RNA genes using MitoAnnotator in MitoFish version 3.72 (Iwasaki et al., 2013; Sato et al., 2018). 
2. Input the gene feature file from MitoAnnotator into Sequence Manipulation Suite 2 (Stothard, 2000) (set to Table 2 - vertebrate mitochondrial DNA) to check that the reading frame is correct for each PCG (no internal stop codons). 
3. Follow the guidelines for submission to GenBank using Bankit (https://www.ncbi.nlm.nih.gov/WebSub/). After submission, download and save the GenBank files in this folder and import them into Geneious (version 2019.1.3).
3. Check the annotated sequences in Geneious to ensure completeness and to manually count overlapping regions and intergenic spaces between PCGs, rRNAs, tRNAs, and non-coding regions.
4. Calculate A+T and G+T content and relative synonymous codon usage (RSCU) of PCGs in DAMBE v. 7.0.35 (Xia, 2001). Base composition skewness formula: AT-skew = [A-T]/[A + T] and GC-skew = [G-C]/[G + C] (Perna and Kocher, 1995). 
5. Make graphs for nucleotide composition and RSCU in R (see scripts in relevant folders).
6. Predict tRNA secondary structure using the generalized vertebrate mitochondrial tRNA settings in ARWEN v. 1.2.3 (Björn Canbäck Bioinformatics) (Laslett and Canbäck, 2008) and the tRNAscanSE webserver v. 2.0 (http://lowelab.ucsc.edu/cgi-bin/tRNAscan-SE2.cgi) (Lowe and Chan, 2016).
7. Characterise control region repetitive regions using the “Tandem Repeat Finder” webserver (https://tandem.bu.edu/trf/trf.html) (Benson, 1999) maintaining default settings.
8. Download the graphic produced by MitoAnnotator and save the files as a png. Edit and enlarge gene names and incude species-specific images, gene numberd and total length in the center.

## STEP 4: Sequence alignments and concatenation

### STEP 4.1: Ingroup and outgroup mitogenome retrieval from GenBank 

1. Save the five newly assembled houndshark mitogenomes in Genbank (full) format in 1a_Data.
2. Compile a list of mitogenome accession numbers from Wang et al. 2022 and Kousteni et al. 2021 as well as four outgroups each from the Lamniform and Orectolobiform orders and save the text file as kousteni-wang_mitogenomes_genbank.list.
3. Use kousteni-wang_mitogenomes_genbank.list in a Batch Entrez(https://www.ncbi.nlm.nih.gov/sites/batchentrez) search to retrieve and download the mitogenome records in Genbank (full) format as a file named kousteni-wang.gb.
4. Merge kousteni-wang.gb and the five newly assembled mitogenomes:
cat *.gb > winn_2022.gb.

### STEP 4.2: Gene region extrations

1. Extract rRNA and CDS DNA sequences from the Genbank file using GBSEQEXTRACTOR v0.04 (https://github.com/linzhi2013/gbseqextractor).
```
gbseqextractor -f winn_2022.gb -prefix winn_2022 -types rRNA -s # output file ´winn_2022.rrna.fasta´
gbseqextractor -f winn_2022.gb -prefix winn_2022 -types CDS -s	# output file `winn_2022.cds.fasta´
```

2. Merge the rRNA and CDS fasta files.
```
cat winn_2022.rrna.fasta winn_2022.cds.fasta > winn_2022.cds-rrna.fasta
	# Using Notepad+++, edit the file ´winn_2022.cds-rrna.fasta´ to standardize gene names.
	# For instance, some GenBank records call 12S rRNA as s-rRNA, 12S ribosomal RNA or rrnS,
	# CO1 as COX1, ND2 as nad2 etc. Standardize all genes to the following code: 
	# 'ATP6' 'ATP8' 'COX1' 'COX2' 'COX3' 'CYTB' 'ND1' 'ND2' 'ND3' 'ND4' 'ND4L' 'ND5' 'ND6' '12SrRNA' '16SrRNA'.
	# Save the edited file as `winn_2022.cds-rrna.std.fasta´.
```

3. Extract individual gene sequences (.fa) from winn_2022.cds-rrna.std.fasta using the custom script maduna2022-gene-extractions.sh
```
#!/bin/bash

GENE=('ATP6' 'ATP8' 'COX1' 'COX2' 'COX3' 'CYTB' 'ND1' 'ND2' 'ND3' 'ND4;' 'ND4L' 'ND5' 'ND6' '12SrRNA' '16SrRNA')

for g in "${GENE[@]}"
	do
		awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) \
		{ printf("%s", $0); } else { printf("\t%s", $0); } }' \
		winn_2022.cds-rrna.std.fasta | grep -F  $g - | tr "\t" "\n" > "${g}".fa
	done

mv 'ND4;.fa' ND4.fa

for f in *.fa
	do
		sed -i 's/;/_/g' $f
	done
```

### STEP 4.3: Multiple Sequence Alignment

1. Create a new folder called 2a_MACSE2. Download the latest version of MACSE and save 'macse_v2.06.jar' at present, in the folder.
2. Save the extracted protein coding (PCG) genes from 1b_GeneExtract, first checking that they are in the correct orientation in Geneious 2019.2.1, as .fasta files in 2a_MACSE2.
3. Align the PCGs using the for loop script maduna2022-13pcgs-msa.sh.
```
datadir=./2a_MACSE2
for i in $datadir/*.fa
do		
	java -jar -Xmx600m macse_v2.06.jar -prog alignSequences -seq "$i" -gc_def 2		
done
	# alignSequences: aligns nucleotide (NT) coding sequences using their amino acid (AA) translations.
	# gc_def: specify the genetic code 2 (The_Vertebrate_Mitochondrial_Code) or change accordingly for your study taxa.
```
4. Align the rRNA genes using the online version of MAFFT (version 7, https://mafft.cbrc.jp/alignment/server/) using the Q-INS-i iterative refinement method, adjusting the direction according to the first sequences confirmed to be in the 5' to 3' direction.
5. Open all alignments from 2a_MACSE in Geneious 2019.2.1 and remove stop codons and ensure each alignment has a length divisible by 3.
6. Change the translation setting to Vertebrate Mitochondrial (Table 2) and check for internal stop codons. 
7. If there are remaining ambiguosly aligned sites, remove them with BMGE version 1.12_1 maintaining default settings through the NGPhylogeny.fr webserver (https://ngphylogeny.fr/).
8. Clean both rRNA alignments with BMGE.
9. Export the edited alignments into the folder 2c_CleanEdit
10. Concatenate the 13 PCGs in Geneious 2019.2.1 and save as 13PCGs_NT (Dataset 1) in fasta, nexus and phyllip format.
11. Concatenate the 13 PCGs and 2 rRNA genes in Geneious and save as 13PCGs_2rRNAs_NT (Dataset 2) in fasta, nexus and phyllip format.
12. Translate 13PCGs_NT and save as 13PCGs_AA (Dataset 3) in fasta, nexus and phyllip format.
13. Make length summaries with the length and alignment locations of each alignment to use for the partition files.

## STEP 5: Substitution saturation and data partitioning schemes 

### STEP 5.1: construct 10 partition nexus files, two for dataset 1, six for dataset 2 and two for dataset 3.

Dataset 1: one partition for the entire alignement, 13 paritions for each PCG.
Dataset 2: one partition for entire alignment, two paritions (rRNA and PCGs), 15 partitions (for each rRNA and PCG), 28 partitions (for each rRNA and for codon position 1 and 2), 28 partitions (for each rRNA and for codon position 1 and 3), 41 partitions (for each rRNA and for each codon position in each PCG).
Dataset 3: one partition for the entire alignment, 13 partitions for each PCG.
***[See "Partitions" to access these files]

### STEP 5.2: test for nucleotide saturation

1. Perform two-tailed tests to examine the degree of nucleotide substitution saturation (Xia et al., 2003) for each gene as well as each codon position of the 13 PCGs, 13PCGs_NT and 13PCGs_rRNAs_NT.
Take into account the proportion of invariant sites as recommended by Xia and Lemey (2009) in DAMBE v.7.2.141 (Xia, 2018). 
2. Visually inspect substitution saturation by plotting the number of transitions (s) and transversions (v) versus divergence.
Divergence is based on genetic distances derived from the Kimura two-parameter (K2P or K80) substitution model (Kimura, 1980). 
The K80 substitution model accommodates transition/transversion rate. 
Bias and the ‘K80 distance’ is expected to increase linearly with divergence time. 

### STEP 5.3: use MODELFINDER v. 1.6.12 (Kalyaanamoorthy et al., 2017) in IQTree v. 2.1.3 (Minh et al., 2020) to determine the best partitioning scheme & corresponding evolutionary models to use in a Maximum Likelihood analysis.

We chose the new model selection procedure (− m MF + MERGE), which additionally implements the FreeRate heterogeneity model inferring the site rates directly from the data instead of being drawn from a gamma distribution (Soubrier et al., 2012). 
The top 30% partition schemes were checked using the relaxed clustering algorithm (− rcluster 30), as described in Lanfear et al. (2014). 
```
#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./3_ML
for i in $datadir/*.nex
	-s "$i"
	-spp $datadir/13PCGs_NT 
	-m MF+MERGE 
	-rcluster 30 
	-AICc
	-nt AUTO
	--safe
done

#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./3_ML
for i in $datadir/*.nex
	-s "$i"
	-spp $datadir/13PCGs_rRNAs_NT 
	-m MF+MERGE 
	-rcluster 30 
	-AICc
	-nt AUTO
	--safe
done

#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./3_ML
for i in $datadir/*.nex
	-s "$i"
	-spp $datadir/13PCGs_AA
	-m MF+MERGE 
	-rcluster 30 
	-AICc
	-nt AUTO
	--safe
done
```
### STEP 5.4: Use MODELFINDER in IQTree to determine the best partitioning scheme & corresponding evolutionary models to use in BI analyses
Use mset mrbayes to restrict the results to models supported by MrBayes.
```
#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./4_BI
for i in $datadir/*.nex
	-s "$i"
	-spp $datadir/13PCGs_NT 
	-m MF+MERGE 
	-rcluster 30 
	-AICc
	-mset mrbayes
	-nt AUTO
	--safe
done

#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./4_BI
for i in $datadir/*.nex
	-s "$i"
	-spp $datadir/13PCGs_rRNAs_NT 
	-m MF+MERGE 
	-rcluster 30 
	-AICc
	-mset mrbayes
	-nt AUTO
	--safe
done

#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./4_BI
for i in $datadir/*.nex
	-s "$i"
	-spp $datadir/13PCGs_AA
	-m MF+MERGE 
	-rcluster 30 
	-AICc
	-mset mrbayes
	-nt AUTO
	--safe
done
```
### STEP 5.5: evaluate the nucleotide substitution values and AICc to select the best partitioning scheme.
The initial partitioning scheme with the lowest AIC and highest concordance factor (straight line graph) should be the most accurate.

## STEP 6: Phylogenetic reconstructions 

### STEP 6.1: use the best-fitting partitioning scheme to construct a ML tree
```
#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./3_ML
	-s $datadir/13PCGs_rRNAs_NT
	-spp $datadir/12_3_NT.nex.best_scheme.nex 
	-alrt 1000 
	-bb 1000
	-nt AUTO
	--safe
done

#!/bin/bash
module load app/IQTREE/1.6.12
iqtree
datadir=./3_ML
	-s $datadir/13PCGs_AA 
	-spp $datadir/genes_AA.nex.best_scheme.nex 
	-alrt 1000 
	-bb 1000
	-nt AUTO
	--safe
done
```

### STEP 6.2: BI tree reconstructio n

### STEP 6.3: constructing the consensus tree

Open the nex.treefile from the ML analysis in Figtree and root it at the outgroup. Save the tree (with bootstrap values) in newick format.
Navigate the the Evolview v3 webpage (https://www.evolgenius.info/evolview/) and make a new project.
Import the newick file.
Adjust size and layout and select bootstrap values.
Import annotations. We used the following:

#### Family groups
```
## style 1	leaf_name	text=mammal,color=darkgreen,textorientation=vertical,linewidth=4,fontsize=16,linestyle=dashed
## style 2-5	leaf_name	bkcolor=#BE4144,text=mammal,textorientation=vertical,linewidth=4,fontsize=16
!grouplabel	style=2
!op	0.8
Chiloscyllium_griseum_NC_017882,Orectolobus_japonicus_NC_022148	bkcolor=darkgrey,text=Orectolobiformes
Carcharhinus_brachyurus_NC_057525,Carcharhinus_leucas_NC_023522	bkcolor=#FCDABE,text=Carcharhinidae
Sphyrna_lewini_NC_022679,Sphyrna_zygaena_NC_025778	bkcolor=#B7F1A5,text=Sphyrnidae
Hemigaleus_microstoma_NC_029400,Hemipristis_elongata_NC_032065	bkcolor=#FF9AA2,text=Hemigaleidae
Hemitriakis_japanica_NC_026774,Galeorhinus_galeus_ON652874	bkcolor=#B2EBF9,text=Triakidae
Proscyllium_habereri_NC_030216	bkcolor=#FFFEC4,text=Proscyllidae
Pseudotriakis_microdon_NC_022735	bkcolor=#F9D2EF,text=Pseudotriakidae
Cephaloscyllium_umbratile_NC_029399,Poroderma_pantherinum_NC_043830	bkcolor=#C6B7F1,text=Scyliorhinidae
Alopias_pelagicus_NC_022822, Lamna_ditropis_NC_024269	bkcolor=lightgrey,text=Lamniformes
```
#### Mustelus genus reproductive mode
```
!grouplabel	style=1
!op	0.8
Mustelus_asterias_ON652873,Triakis_megalopterus_ON075075	text=aplacental-spotted,color=#0033CC,textorientation=horizontal,linewidth=4,fontsize=14,linestyle=dashed
Mustelus_griseus_NC_023527,Mustelus_mustelus_NC_039629	text=placental-non-spotted,color=#3399FF,textorientation=horizontal,linewidth=4,fontsize=14 
```
#### Bootstrap style
```
!bootstrapValueStyle	show=2,style=numeric,color=darkred,place=2,size=12
```
#### Marking our sequences
```
Mustelus_asterias_ON652873	star,red	
Mustelus_palumbes_ON075076	star,red
Triakis_megalopterus_ON075075	star,red
Mustelus_mosis_ON075077		star,red
Galeorhinus_galeus_ON652874	star,red
```
