# Bioinformatics pipeline of Winn _et al_. (2023) - Triakid Mitophylogenomics

## Background 

Complex evolutionary patterns in the mitochondrial genome (mitogenome) of the most species-rich order, the Carcharhiniforms (groundsharks) has yielded challenges in phylogenomic reconstruction of families and genera belonging to it, particularly in the family Triakidae (houndsharks), where there are arguments for both monophyly and paraphyly. We hypothesized that opposing resolutions are a product of the _a priori_ partitioning scheme selected. Accordingly, we employed an extensive statistical framework to select our partitioning scheme for inference of the mitochondrial phylogenomic relationships within Carcharhiniforms and used the multispecies coalescent model to account for the influence of gene tree discordance on species tree inference. We included five new houndshark mitogenomes to increase resolution of Triakidae and uncovered a 714 bp-duplication in the assembly of _Galeorhinus galeus_. Phylogenetic reconstruction confirmed monophyly within Triakidae and the existence of two clades of the expanded _Mustelus_ genus, alluding to the evolutionary reversal of reproductive mode from placental to aplacental. 

Here we present, for the first time, the Ion Torrent® next-generation sequencing (NGS) data and the complete mitochondrial genomes (mitogenomes) of five houndsharks (Chondrichthyes: Triakidae), which include _Galeorhinus galeus_ (17,488 bp; GenBank accession number ON652874), _Mustelus asterias_ (16,708; ON652873), _Mustelus mosis_ (16,755; ON075077), _Mustelus palumbes_ (16,708; ON075076), and _Triakis megalopterus_ (16746; ON075075). We also present the code used to assemble and annotate their mitogenomes, prepare alignments, partition our datasets, assign models of evolution, infer phylogenies based on traditional site homogenous concatenation approaches as well as under the multispecies coalescent model (MSCM) and site heterogenous models, and generate statistical data for comparison of different topological outcomes. The data and code presented in this paper can be used by other researchers to delve deeper into the phylogenetic relationships of Carcharhiniformes (groundsharks) and the diversification of triakid species as mitogenomes accumulate in public repositories.

## STEP 1: Quality control of Ion GeneStudio™ S5 sequencing data

1. Check sequence quality in FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Trim adaptors and poor-quality bases (phred score below 20) and remove reads shorter than 25 base pairs (bp) in Torrent Suite Version 5.16.[^1]

## STEP 2: Mitogenome assembly

### STEP 2.1: Reference-based assembly

1. Align raw reads to the _Mustelus mustelus_ mitogenome (NC_039629.1) using the Geneious read mapper with medium sensitivity settings and five iterations in Geneious Prime (version 2019.1.3).
[^1]: Trimmed and reference-filtered sequencing reads in BAM format can be found in the "1_Raw Ion Torrent® NGS data" folder.

### STEP 2.2: Hybrid assembly

1. Extract reads that mapped to the reference genome in bam format.
2. Feed the reads into a _de novo_ pipeline in SPAdes version 3.15 with the input set for unpaired Ion Torrent reads with 8 threads, kmers 21,33,55,77,99,127, the careful option to reduce the number of mismatches and short indels and all other parameters left as default.
EXAMPLE - _Mustelus palumbes_

```bash

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

### STEP 2.3: _de novo_ assembly

1. Directly map raw reads (bam format) _de novo_.
EXAMPLE - _Mustelus palumbes_

```bash
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
2. Check each alignment for discrepancies on  Geneious and edit maually to obtain the final genome sequence.
3. If there is a significant discrepancy between the three alignments further investigation is warranted. We Sanger sequenced a region of our _Galeorhinus galeus_ mitogenome that we could not find a consensus on when comparing our three assemblies. We compared the Sanger sequence fragment to our three assemblies in Geneious and found that it matched a section of the _de novo_ assembly that the reference assembly failed to detect. The full details on this step can be found in the Supplementary material of Winn et al. (2023).

## STEP 3: Mitogenome annotation

1. Annotate protein coding genes (PCGs), ribosomal (r)RNA and transfer (t)RNA genes using MitoAnnotator in MitoFish version 3.72 (Iwasaki et al., 2013; Sato et al., 2018). 
2. Input the gene feature file from MitoAnnotator into Sequence Manipulation Suite 2 (Stothard, 2000) (set to Table 2 - vertebrate mitochondrial DNA) to check that the reading frame is correct for each PCG (no internal stop codons). 
3. Follow the guidelines for submission to GenBank using Bankit (https://www.ncbi.nlm.nih.gov/WebSub/). After submission, download and save the GenBank files in this folder and import them into Geneious (version 2019.1.3).
3. Check the annotated sequences in Geneious to ensure completeness and manually count overlapping regions and intergenic spaces between PCGs, rRNAs, tRNAs, and non-coding regions.
4. Calculate A+T and G+T content and relative synonymous codon usage (RSCU) of PCGs in DAMBE v. 7.0.35 (Xia, 2001). Base composition skewness formula: AT-skew = [A-T]/[A + T] and GC-skew = [G-C]/[G + C] (Perna and Kocher, 1995). 
5. Make graphs for nucleotide composition and RSCU in R. [^2]
[^2]: winn2023-nucleotide-composition-plots and winn2023-codon-usage-plots can be found in the "5_Scripts" folder.
8. Predict tRNA secondary structure using the generalized vertebrate mitochondrial tRNA settings in ARWEN v. 1.2.3 (Björn Canbäck Bioinformatics) (Laslett and Canbäck, 2008) and the tRNAscanSE webserver v. 2.0 (http://lowelab.ucsc.edu/cgi-bin/tRNAscan-SE2.cgi) (Lowe and Chan, 2016). We also used the MXfold2 webserver (http://ws.sato-lab.org/mxfold2/) to verify the secondary structure of the duplicate tRNA-Thr in _Galeorhinus galeus_.
9. Characterise control region repetitive regions using the “Tandem Repeat Finder” webserver (https://tandem.bu.edu/trf/trf.html) (Benson, 1999) maintaining default settings.
10. Download the graphic produced by MitoAnnotator and save the files as a png. Edit and enlarge gene names and incude species-specific images, gene numbers and total length in the center.[^3]
[^3]: Fully annotated GenBank files of each of our newly assembled mitogenomes can be found on GenBank with accession numbers ON075075, ON075076, ON075077, ON652873, and ON652874.

## STEP 4: Sequence alignment and concatenation

### STEP 4.1: Ingroup and outgroup mitogenome retrieval from GenBank 

1. Save the five newly assembled houndshark mitogenomes in Genbank (full) format.
2. Compile a list of mitogenome accession numbers from Wang et al. 2022 and Kousteni et al. 2021 as well as four outgroups each from the Lamniform and Orectolobiform orders and save the text file as kousteni-wang_mitogenomes_genbank.list.
3. Use kousteni-wang_mitogenomes_genbank.list in a Batch Entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez) search to retrieve and download the mitogenome records in Genbank (full) format as a file named kousteni-wang.gb.
4. Merge kousteni-wang.gb and the five newly assembled mitogenomes.[^4]
```
cat *.gb > winn_2023.gb.
```
[^4]: GenBank files and accession number lists for our dataset are in the folder "2_Galeomorphii mitogenome sequences". 

### STEP 4.2: Gene region extractions

1. Extract rRNA and CDS DNA sequences from the Genbank file using GBSEQEXTRACTOR v0.04 (https://github.com/linzhi2013/gbseqextractor).
```
gbseqextractor -f winn_2023.gb -prefix winn_2023 -types rRNA -s # output file ´winn_2023.rrna.fasta´
gbseqextractor -f winn_2023.gb -prefix winn_2023 -types CDS -s	# output file `winn_2023.cds.fasta´
```
2. Merge the rRNA and CDS fasta files.
```
cat winn_2023.rrna.fasta winn_2023.cds.fasta > winn_2023.cds-rrna.fasta
```
3. Using Notes, edit the file ´winn_2023.cds-rrna.fasta´ to standardize gene names.
For instance, some GenBank records denote 12S rRNA as s-rRNA, 12S ribosomal RNA or rrnS, CO1 as COX1, ND2 as nad2 etc. 
Standardize all genes to the following code:  'ATP6' 'ATP8' 'COX1' 'COX2' 'COX3' 'CYTB' 'ND1' 'ND2' 'ND3' 'ND4' 'ND4L' 'ND5' 'ND6' '12SrRNA' '16SrRNA'.
4. Save the edited file as `winn_2023.cds-rrna.std.fasta´.
5. Extract individual gene sequences (.fa) from winn_2023.cds-rrna.std.fasta using the custom script maduna2022-gene-extractions.sh:

```bash
#!/bin/bash
GENE=('ATP6' 'ATP8' 'COX1' 'COX2' 'COX3' 'CYTB' 'ND1' 'ND2' 'ND3' 'ND4;' 'ND4L' 'ND5' 'ND6' '12SrRNA' '16SrRNA')
for g in "${GENE[@]}"
	do
		awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) \
		{ printf("%s", $0); } else { printf("\t%s", $0); } }' \
		winn_2023.cds-rrna.std.fasta | grep -F  $g - | tr "\t" "\n" > "${g}".fa
	done
mv 'ND4;.fa' ND4.fa
for f in *.fa
	do
		sed -i 's/;/_/g' $f
	done
```

### STEP 4.3: Multiple Sequence Alignment

1. Create a new folder called 2a_MACSE2. Download the latest version of MACSE and save 'macse_v2.06.jar' at present, in the folder.
2. Check that the extracted protein coding genes  (PCGs) from 1b_GeneExtract are in the correct orientation in Geneious 2019.2.1 and save them as .fasta files in 2a_MACSE2.
3. Align the PCGs using the for loop script maduna2022-13pcgs-msa.sh.

```bash
#!/bin/bash
datadir=./2a_MACSE2
for i in $datadir/*.fa
do		
	java -jar -Xmx600m macse_v2.06.jar -prog alignSequences -seq "$i" -gc_def 2		
done
	# alignSequences: aligns nucleotide (NT) coding sequences using their amino acid (AA) translations.
	# gc_def: specify the genetic code 2 (The_Vertebrate_Mitochondrial_Code) or change accordingly for your study taxa.
```
4. Align the rRNA genes using the online version of MAFFT (version 7, https://mafft.cbrc.jp/alignment/server/) using the Q-INS-i iterative refinement method, adjusting the direction according to the first sequences confirmed to be in the 5' to 3' direction.
5. Open all alignments from 2a_MACSE in Geneious 2019.2.1, change the translation settings to Vertebrate Mitochondrial (Table 2), remove stop codons and ensure each alignment has a length divisible by 3.
6. If there are remaining ambiguosly aligned sites, remove them with BMGE version 1.12_1 maintaining default settings through the NGPhylogeny.fr webserver (https://ngphylogeny.fr/). Clean both the rRNA alignments with BMGE.
7. Export the edited alignments into the folder 2c_CleanEdit
8. Concatenate the 13 PCGs in Geneious 2019.2.1 and save as 13PCGs_NT (Dataset 1) in fasta, nexus and phyllip format.
9. Concatenate the 13 PCGs and 2 rRNA genes in Geneious and save as 13PCGs_2rRNAs_NT (Dataset 2) in fasta, nexus and phyllip format.
10. Translate 13PCGs_NT and save as 13PCGs_AA (Dataset 3) in fasta, nexus and phyllip format.
11. Make and save length summaries with the length and alignment locations of each alignment to use for the partition files.[^5]
[^5]: Concatenated datasets can be found in the folder "3_Multiple sequence alignments".

## STEP 5: Substitution saturation and _a priori_ data partitioning

### STEP 5.1: Test for nucleotide saturation.

1. Perform two-tailed tests to examine the degree of nucleotide substitution saturation (Xia et al., 2003) for each gene and each codon position of the 13PCGs as well as on the entire 13_PCGs_NT and 13PCGs_rRNAs_NT datasets, taking into account the proportion of invariant sites as recommended by Xia and Lemey (2009), in DAMBE v.7.2.141 (Xia, 2018).
3. Visually inspect substitution saturation by plotting the number of transitions (s) and transversions (v) versus divergence.
Divergence is based on genetic distances derived from the Kimura two-parameter (K2P or K80) substitution model (Kimura, 1980). 
The K80 substitution model accommodates transition/transversion rate. 
Bias and the ‘K80 distance’ is expected to increase linearly with divergence time.
See: https://www.researchgate.net/publication/232048505_Assesing_substitution_saturation_with_DAMBE for a more detailed tutorial.

### STEP 5.2: Construct partition nexus files.

Dataset 1 & 3: one partition for the entire alignement, 13 paritions for each PCG.
Dataset 2: one partition for the entire alignment, five paritions (codon: pos1 + pos2 + pos3 + 2 rRNAs ), four partitions (codon: pos1_pos2 + pos3 + 2 rRNAs), four partitions (codon: pos1 + pos3 + 2 rRNAs), 15 partitions for each PCG and 2 rRNAs, 41 partitions (geneXcodon: 13 PCGs pos1 + 13 pos2 + 13 pos2 + 2 rRNAs), 28 partitions (geneXCodon: 13 PCGs pos1_pos2 + 13 PCGs pos3 + 2 rRNAs), 28 partitions (geneXcodon: 13 PCGs pos1 + 13 PCGs pos3 + 2 rRNAs).

## STEP 6: Mitophylogenomic reconstruction[^5][^6]

[^6]: Partition files are available in "4_Partition Files".

### STEP 6.1: Use MODELFINDER v. 1.6.12 (Kalyaanamoorthy et al., 2017) in IQTree v. 2.1.3 (Minh et al., 2020) to determine the best partitioning scheme and corresponding evolutionary models to use in a Maximum Likelihood analysis.

Apply the new model selection procedure (-m MF+MERGE) which additionally implements the FreeRate heterogeneity model, inferring the site rates directly from the data instead of being drawn from a gamma distribution (-cmax 20; Soubrier et al., 2012).
The top 30% partition schemes were checked using the relaxed clustering algorithm (− rcluster 30), as described in Lanfear et al. (2014). 
The maximum number of rate categories, cmax, is set to 20 since it is likely to observe more rate variations for alignments with many sequences. 
3 CPU cores, T, are used to decrease computational burden.

1. 13 PCGs_NT

```bash
#!/bin/bash
module load app/IQTREE/2.1.3
cd *./13PCGs
iqtree2 -s 13PCGs_NT.fasta -st DNA -p PS1/PS1.txt -pre PS1/PS1_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_NT.fasta -st DNA -p PS5/PS5.txt -pre PS5/PS5_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
```
2. 13PCGs_2rRNAs_NT

```bash
#!/bin/bash
module load app/IQTREE/2.1.3
cd *./13PCGs_2rRNAs
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS1/PS1.txt -pre PS1/PS1_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS2/PS2.txt -pre PS2/PS2_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS3/PS3.txt -pre PS3/PS3_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS4/PS4.txt -pre PS4/PS4_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS5/PS5.txt -pre PS5/PS5_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS6/PS6.txt -pre PS6/PS6_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS7/PS7.txt -pre PS7/PS7_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS8/PS8.txt -pre PS8/PS8_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
```
3. 13PCGs_AA

```bash
#!/bin/bash
module load app/IQTREE/2.1.3
cd *./13PCGs
iqtree2 -s 13PCGs_AA.fasta -st AA -p PS1/PS1AA.txt -pre PS1/PS1AA_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
iqtree2 -s 13PCGs_AA.fasta -st AA -p PS5/PS5AA.txt -pre PS5/PS5AA_run01_mf -m MF+MERGE -AICc -rcluster 30 -T 3 -cmax 20
```

Repeat the above using the Bayesian Information Critereon (BIC)

### STEP 6.2: Use MODELFINDER in IQTree to determine the best partitioning scheme & corresponding evolutionary models to use in BI analyses.

Apply secondary model selection for the best-fitting partitioning identified by ModelFinder in Step 6.1 under the FreeRate heterogeneity model to select the next best model for Bayesian inference.
Rerun ModelFinder with options: -m TESTONLY -mset mrbayes to restrict the results to models supported by MrBayes.

1. 13PCGs_NT
```bash
#!/bin/bash
module load app/IQTREE/2.1.3
cd *./13PCGs
iqtree2 -s 13PCGs_NT.fasta -st DNA -p PS1/PS1_run01_mf.best_scheme.nex -pre PS1/PS1_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_NT.fasta -st DNA -p PS5/PS5_run01_mf.best_scheme.nex -pre PS5/PS5_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
```
2. 13PCGs_2rRNAs

```bash
#!/bin/bash
module load app/IQTREE/2.1.3
cd *./13PCGs_2rRNAs
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS1/PS1_run01_mf.best_scheme.nex -pre PS1/PS1_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS2/PS2_run01_mf.best_scheme.nex -pre PS2/PS2_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS3/PS3_run01_mf.best_scheme.nex -pre PS3/PS3_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS4/PS4_run01_mf.best_scheme.nex -pre PS4/PS4_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS5/PS5_run01_mf.best_scheme.nex -pre PS5/PS5_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS6/PS6_run01_mf.best_scheme.nex -pre PS6/PS6_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS7/PS7_run01_mf.best_scheme.nex -pre PS7/PS7_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS8/PS8_run01_mf.best_scheme.nex -pre PS8/PS8_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
```
3. 13PCGs_AA

```bash
#!/bin/bash
module load app/IQTREE/2.1.3
cd *./13PCGs
iqtree2 -s 13PCGs_NT.fasta -st AA -p PS1/PS1AA_run01_mf.best_scheme.nex -pre PS1/PS1AA_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
iqtree2 -s 13PCGs_NT.fasta -st AA -p PS5/PS5AA_run01_mf.best_scheme.nex -pre PS5/PS5AA_run01_mf -m TESTONLY -AICc -T AUTO -mset mrbayes
```

### STEP 6.3: Evaluate likelihood statistics to select the best partitioning scheme.

These values can be found in the IQTREE files generated by ModelFinder.
The partitioning scheme with the lowest AICc/BIC and highest concordance factor (straight line graph) should be the most accurate.

### STEP 6.4: Construct ML trees for each partitioning scheme.

Use the substitution models indicated in best_model.nex files for each partitioning scheme to construct Maximum Likelihood phylogenies.
Use the Nearest Neighbor Interchange (NNI) approach to search for tree topology.
Compute branch supports with 1000 replicates of the Shimodaira-Hasegawa approximate likelihood-ratio test (SH-aLRT; Anisimova and Gascuel, 2006) and the ultrafast bootstrapping (UFBoot2) approach (Hoang et al., 2018).

1. 13PCGs_NT

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs
iqtree2 -s 13PCGs_NT.fasta -st DNA -p PS1/PS1_run01_mf.best_model.nex -pre PS1/PS1_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_NT.fasta -st DNA -p PS5/PS5_run01_mf.best_model.nex -pre PS5/PS5_run02_ml -T 3 -B 1000 -alrt 1000
```
2. 13PCGs_2rRNAs_NT

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs_2rRNAs
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS1/PS1_run01_mf.best_model.nex -pre PS1/PS1_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS2/PS2_run01_mf.best_model.nex -pre PS2/PS2_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS3/PS3_run01_mf.best_model.nex -pre PS3/PS3_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS4/PS4_run01_mf.best_model.nex -pre PS4/PS4_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS5/PS5_run01_mf.best_model.nex -pre PS5/PS5_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS6/PS6_run01_mf.best_model.nex -pre PS6/PS6_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS7/PS7_run01_mf.best_model.nex -pre PS7/PS7_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -st DNA -p PS8/PS8_run01_mf.best_model.nex -pre PS8/PS8_run02_ml -T 3 -B 1000 -alrt 1000
```
3. 13PCGs_AA

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs
iqtree2 -s 13PCGs_AA.fasta -st AA -p PS1/PS1AA_run01_mf.best_model.nex -pre PS1/PS1AA_run02_ml -T 3 -B 1000 -alrt 1000
iqtree2 -s 13PCGs_AA.fasta -st AA -p PS5/PS5AA_run01_mf.best_model.nex -pre PS5/PS5AA_run02_ml -T 3 -B 1000 -alrt 1000
```

### STEP 6.5: Perform Bayesian Inference analysis using the Cyberinfrastructure for Phylogenetic Research (CIPRES) Science Gateway portal v. 3.3 (www.phylo.org) at the San Diego Supercomputer Center (Miller et al., 2010).[^7]

1. Run a pair of independent searches for 5 million generations, with trees saved every 1,000 generations and the first 2,500 sampled trees of each search discarded as burn-in. 
2. Screen the model parameter summary statistics Estimated Sample Size (ESS) and Potential Scale Reduction Factor (PSRF), where convergence occurred at ESS >200, and PSRF ~1.0. 

```bash
#!/bin/bash
for g in *.nex
	do
		mb -i $g
	done
```
[^7]: Scripts to plot MrBayes log files are stored as kmisc.R, mb_plots.R in "5_Scripts" - https://rdrr.io/github/kmiddleton/kmmisc/man/plot_mrb.html.

### STEP 6.6: Compute confordance factors

Investigate topological conflict around each branch of the species tree by calculating gene and site concordance factors in IQ-Tree.
Infer concatenation-based species trees with 1000 ultrafast bootstraps and an edge-linked partition model.

1. 13PCGs_NT

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs

# Calculate gene concordance factors (gCF).
iqtree2 -s 13PCGs_NT.fasta -p PS1/PS1_run01_mf.best_scheme.nex --prefix PS1/PS1_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_NT.fasta -p PS5/PS5_run01_mf.best_scheme.nex --prefix PS5/PS5_run03_concat.condonpart.MF -B 1000 -T 3

# Calculate site concordance factors (sCF) and infer the locus trees.
iqtree2 -s 13PCGs_NT.fasta -S PS1/PS1_run01_mf.best_scheme.nex --prefix PS1/PS1_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_NT.fasta -S PS5/PS5_run01_mf.best_scheme.nex --prefix PS5/PS5_run03_loci.condonpart.MF -T 3

# Compute concordance factors.
iqtree2 -t PS1/PS1_run03_concat.condonpart.MF.treefile --gcf PS1/PS1_run03_loci.condonpart.MF.treefile -s 13PCGs_NT.fasta --scf 100 -seed 471990 --prefix PS1/PS1_run03_concord
iqtree2 -t PS5/PS5_run03_concat.condonpart.MF.treefile --gcf PS5/PS5_run03_loci.condonpart.MF.treefile -s 13PCGs_NT.fasta --scf 100 -seed 471990 --prefix PS5/PS5_run03_concord
```

2. 13PCGs_2rRNAs_NT

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs_2rRNAs

# Calculate gene concordance factors (gCF).
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS1/PS1_run01_mf.best_scheme.nex --prefix PS1/PS1_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS2/PS2_run01_mf.best_scheme.nex --prefix PS2/PS2_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS3/PS3_run01_mf.best_scheme.nex --prefix PS3/PS3_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS4/PS4_run01_mf.best_scheme.nex --prefix PS4/PS4_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS5/PS5_run01_mf.best_scheme.nex --prefix PS5/PS5_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS6/PS6_run01_mf.best_scheme.nex --prefix PS6/PS6_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS7/PS7_run01_mf.best_scheme.nex --prefix PS7/PS7_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -p PS8/PS8_run01_mf.best_scheme.nex --prefix PS8/PS8_run03_concat.condonpart.MF -B 1000 -T 3

# Calculate site concordance factors (sCF) and infer the locus trees.
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS1/PS1_run01_mf.best_scheme.nex --prefix PS1/PS1_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS2/PS2_run01_mf.best_scheme.nex --prefix PS2/PS2_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS3/PS3_run01_mf.best_scheme.nex --prefix PS3/PS3_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS4/PS4_run01_mf.best_scheme.nex --prefix PS4/PS4_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS5/PS5_run01_mf.best_scheme.nex --prefix PS5/PS5_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS6/PS6_run01_mf.best_scheme.nex --prefix PS6/PS6_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS7/PS7_run01_mf.best_scheme.nex --prefix PS7/PS7_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -S PS8/PS8_run01_mf.best_scheme.nex --prefix PS8/PS8_run03_loci.condonpart.MF -T 3

# Compute concordance factors.
iqtree2 -t PS1/PS1_run03_concat.condonpart.MF.treefile --gcf PS1/PS1_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS1/PS1_run03_concord
iqtree2 -t PS2/PS2_run03_concat.condonpart.MF.treefile --gcf PS2/PS2_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS2/PS2_run03_concord
iqtree2 -t PS3/PS3_run03_concat.condonpart.MF.treefile --gcf PS3/PS3_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS3/PS3_run03_concord
iqtree2 -t PS4/PS4_run03_concat.condonpart.MF.treefile --gcf PS4/PS4_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS4/PS4_run03_concord
iqtree2 -t PS5/PS5_run03_concat.condonpart.MF.treefile --gcf PS5/PS5_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS5/PS5_run03_concord
iqtree2 -t PS6/PS6_run03_concat.condonpart.MF.treefile --gcf PS6/PS6_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS6/PS6_run03_concord
iqtree2 -t PS7/PS7_run03_concat.condonpart.MF.treefile --gcf PS7/PS7_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS7/PS7_run03_concord
iqtree2 -t PS8/PS8_run03_concat.condonpart.MF.treefile --gcf PS8/PS8_run03_loci.condonpart.MF.treefile -s 13PCGs_2rRNAs_NT.fasta --scf 100 -seed 471990 --prefix PS8/PS8_run03_concord
```

3. 13PCGs_AA

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs

# Calculate gene concordance factors (gCF).
iqtree2 -s 13PCGs_AA.fasta -p PS1/PS1AA_run01_mf.best_scheme.nex --prefix PS1/PS1AA_run03_concat.condonpart.MF -B 1000 -T 3
iqtree2 -s 13PCGs_AA.fasta -p PS5/PS5AA_run01_mf.best_scheme.nex --prefix PS5/PS5AA_run03_concat.condonpart.MF -B 1000 -T 3

# Calculate site concordance factors (sCF) and infer the locus trees.
iqtree2 -s 13PCGs_AA.fasta -S PS1/PS1AA_run01_mf.best_scheme.nex --prefix PS1/PS1AA_run03_loci.condonpart.MF -T 3
iqtree2 -s 13PCGs_AA.fasta -S PS5/PS5AA_run01_mf.best_scheme.nex --prefix PS5/PS5AA_run03_loci.condonpart.MF -T 3

# Compute concordance factors.
iqtree2 -t PS1/PS1AA_run03_concat.condonpart.MF.treefile --gcf PS1/PS1AA_run03_loci.condonpart.MF.treefile -s 13PCGs_AA.fasta --scf 100 -seed 471990 --prefix PS1/PS1AA_run03_concord
iqtree2 -t PS5/PS5AA_run03_concat.condonpart.MF.treefile --gcf PS5/PS5AA_run03_loci.condonpart.MF.treefile -s 13PCGs_AA.fasta --scf 100 -seed 471990 --prefix PS5/PS5AA_run03_concord
```
### STEP 6.7: Run the approximately unbiased (AU) tree topology test.

Compare the trees from the eight runs to determine significant diﬀerences with the approximately unbiased (AU) tree topology test (Shimodaira, 2002) also implemented in IQ-Tree.
1. Save the final eight treefiles as a list in Newick format.
2. Compute the log likelihood of the set of trees (-z).
3. Set the number of search iterations is set to 0 (model parameters are quickly estimated from an initial parsimony tree).
4. Run tree topology tests using the RELL approximation (Kishino et al., 1990). zb specifies the number of RELL replicates.
5. Perform weighted KH and weighted SH tests (-zw).
6. Conduct an approximately unbiased (AU) test (Shimodaira, 2002) (-au).

```bash
#!/bin/bash
module load app/IQTREE/1.6.12
cd ./13PCGs_2rRNAs
iqtree2 -s 13PCGs_2rRNAs_NT.fasta -z TopoTest_PS1-PS8.treesls --prefix TopoTest_run01 -n 0 -zb 10000 -zw -au -nt AUTO
```
A tree is rejected if its p-value < 0.05 (marked with a - sign) for KH, SH and AU tests.
_bp-RELL_ and _c-ELW_ return posterior weights which are not p-values. The weights sum up to 1 across the trees tested.

### STEP 6.8: Test the hypothesis of equal frequencies.

1. Conduct a χ2-test to determine whether the frequency of gene trees (gCF) and sites (sCF) supporting the two alternative topologies differ significantly as implemented in Lanfear’s R script (Minh et al., 2020) in R v.4.1.2 (R Core Team, 2021).
2. Construct plots comparing gCF and sCF values with UFBoot2 values using Lanfear’s R script. [^8]
[^8]: winn2023-concordance-factors-plots script is in the folder "5_Scripts". 

### STEP 6.9: Visualise and analyse the trees.

1. Open all nex.treefile from the ML and BI analyses in FigTree and root at the outgroup.
2. Align nodes and view in increasing order.
3. Highlight the Triakidae family branches and nodes.
4. Add bootstrap, concordance factors and posterior probability values.
5. Save as png files and visualise.
***[Nodes with UFBoot2 ≥ 95, PP ≥ 95 and SH-aLRT ≥ 80 were considered well supported (Minh et al., 2020b).]
***[Nodes with sCF values below 33% and gCF values that are lower than sCF values or near zero require further attention].

### STEP 6.10: Contruct a consensus tree

1. Save the best supported trees above (with bootstrap values) in newick format.
2. Navigate to the Evolview v3 webpage (https://www.evolgenius.info/evolview/) and make a new project.
3. Import the newick file.
4. Adjust size and layout and select bootstrap values.
5. Import annotations.[^9]
[^9]: Annotations are provided as winn2023-evolview-tree-annotations in "5_Scripts".

## STEP 7: Multispecies coalescent modelling

Estimates the effects of gene-tree conflict on species-tree inference.

### STEP 7.1: ASTRAL

1. Estimate individual gene trees for the 13 PCGs and 2 rRNAs based on the ML criterion in IQ-Tree used the cleaned and edited single gene alignments created in Step 4.3.
Use a greedy model selection strategy (-m MFP)and the NNI approach to search for tree topology and compute branch supports with 1000 bootstrapped replicates of the UFBoot2 approach (Hoang et al., 2018).
```
for gene in *.fas
	do
		iqtree2 -s $gene -st DNA -m MFP -AICc -nt AUTO -B 1000 
	done
```
2. Create a combined file containing all the gene trees.`
```
cat *.treefile > elasmo.15gene.tre
````
3. Convert the treefile to Newick format in FigTree.
4. Run the combined gene tree Newick file through ASTRAL v.5.6.3 (Zhang et al., 2018) using the default options.
```
java -jar astral.5.7.8.jar -i elasmo-mitophy-15G.tre -o elasmo-mitophy-15G-ASTRAL.tre
````
5. Collapse branches with low support using newick utilities and then run ASTRAL.
```
nw_ed  elasmo-mitophy-15G.tre 'i & b<=10' o > elasmo-mitophy-15G-BS10.tre
java -jar astral.5.7.8.jar -i elasmo-mitophy-15G-BS10.tre -o elasmo-mitophy-15G-BS10-ASTRAL.tre
```
6. Annotate the branches of the species tree generated in step 2a.
```
# Full annotation
java -jar astral.5.7.8.jar -q elasmo-mitophy-15G-ASTRAL.tre -i elasmo-mitophy-15G.tre -t 2 -o elasmo-mitophy-15G-scored-t2.tre

# Generate the posterior probabilities for branches of the main topology and the two alternatives which add up to three for each branch.
java -jar astral.5.7.8.jar -q elasmo-mitophy-15G-ASTRAL.tre -i elasmo-mitophy-15G.tre -t 4 -o elasmo-mitophy-15G-scored-t4.tre

# Show quartet support for the main topology and the two alternative topologies.
java -jar astral.5.7.8.jar -q elasmo-mitophy-15G-ASTRAL.tre -i elasmo-mitophy-15G.tre -t 8 -o elasmo-mitophy-15G-scored-t8.tre

# Test for polytomies with the method of Sayyari and Mirarab (2018), which is based on a Chi-Square test among quartet frequencies for nodes, implemented with the -t 10 command. 
java -jar astral.5.7.8.jar -q elasmo-mitophy-15G-ASTRAL.tre -i elasmo-mitophy-15G.tre -t `10 -o elasmo-mitophy-15G-scored-t10.tre
```

### STEP 7.2: SVDQuartets

1. Create a nexus file using Dataset 2: Galeomorphii_13PCGs_2rRNAs_NT.[^5] Use the gene partitions from PS05.[^6]
2. Open the nexus file in PAUP* v4.0a 169 (Swofford, 2003).
3. Implement the multispecies coalescent tree model with random quartet sampling of 100,000 replicates and 1,000 bootstrap replicates (https://phylosolutions.com/tutorials/svdq-qage/svdq-qage-tutorial.html).

## STEP 8: Site heterogenous models

### STEP 8.1: Test for compositional heterogeneity among lineages with AliGROOVE 

Use AliGROOVE to test for compositional heterogeneity among lineages in the nucleotide (Data_14_Galeomorphii_13PCGs_NT.fasta and Data_17_Galeomorphii_13PCGs_2rRNAs_NT.fasta) and amino acid datasets (Data_20_Galeomorphii_13PCGs_AA.fasta). For nucleotide alignments, run the software with the -N option to treat indels as ambiguous characters and tested it without invoking -N so that indels were treated as a fifth character trait.
```
aligroove = ./bin/AliGROOVE

# Without invoking the -N option: indels are treated as 5th character
$aligroove/AliGROOVE_v.1.08.pl -i Data_14_Galeomorphii_13PCGs_NT.fasta
$aligroove/AliGROOVE_v.1.08.pl -i Data_17_Galeomorphii_13PCGs_2rRNAs_NT.fasta
$aligroove/AliGROOVE_v.1.08.pl -i Data_20_Galeomorphii_13PCGs_AA.fasta

# With -N option: indels are treated as ambiguous characters (only possible for nucleotide alignments)
$aligroove/AliGROOVE_v.1.08.pl -N -i Data_14_Galeomorphii_13PCGs_NT.fasta
$aligroove/AliGROOVE_v.1.08.pl -N -i Data_17_Galeomorphii_13PCGs_2rRNAs_NT.fasta
```
### STEP 8.2: Test for compositional heterogeneity among lineages with BaCoCa

Run BaCoCa with default settings, invoking the -r option of the program to generate heat maps in combination with hierarchical clustering.
```
Bacoca = ./BaCoCa_v1.1_Skript_Perl
perl $bacoca/BaCoCa.v1.109.r.pl -i 13PCGs_2rRNAs_NT.fasta -P 13PCGs_2rRNAs_NT.part.txt -r
```
### STEP 8.3: Test different models for phylogenetic reconstruction with PHYLOBAYES_MPI 

Conduct BI analyses using the pb_mpi program of the PHYLOBAYES_MPI package. Run two independent Markov chain Monte Carlo (MCMC) chains per model (GTR + G4, CAT + Poisson + G4, CAT + GTR + G4, and C60 + GTR + G4 for amino acids only) for each supermatrix (Data_15_Galeomorphii_13PCGs_NT.phy, Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy, Data_21_Galeomorphii_13PCGs_AA.phy) for as long as time allows for your project. We ran our chains for several days before checking convergence. Convergence of independent runs can be assessed using maximum difference (maxdiff) in all bipartitions and effective sample size (ESS) following the guidelines specified in the PHYLOBAYES_MPI manual (https://github.com/bayesiancook/phylobayes/blob/master/pbManual4.1.pdf). The bpcomp program is used to generate outputs of the largest (maxdiff) and mean (meandiff) discrepancies observed across all bipartitions, while the tracecomp program is used to estimate the ESS and discrepancy for the summary variables of the models. For both programs, we used a conservative burn-in of 20% of the length of MCMC chain.
```
phylobayes = ./data
mod_sel = ./scripts 
```
1. The GTR model
```
# Dataset 1: 13PCGs_NT
mpirun -np 8 $phylobayes/pb_mpi -ncat 1 -gtr -d Data_15_Galeomorphii_13PCGs_NT.phy GTR_13PCGs_NT_chain1 &
mpirun -np 8 $phylobayes/pb_mpi -ncat 1 -gtr -d Data_15_Galeomorphii_13PCGs_NT.phy 
GTR_13PCGs_NT_chain2 &
$phylobayes/tracecomp GTR_13PCGs_NT_chain?.trace -o GTR_13PCGs_NT_tracecomp
$phylobayes/bpcomp -x 813 10 GTR_13PCGs_NT_chain1 GTR_13PCGs_NT_chain2 -o GTR_13PCGs_NT_bpcomp

# Dataset 2: 13PCGs_rRNAs_NT
mpirun -np 8 $phylobayes/pb_mpi -ncat 1 -gtr -d Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy GTR_13PCGs_2rRNAs_NT_chain1 &
mpirun -np 8 $phylobayes/pb_mpi -ncat 1 -gtr -d Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy GTR_13PCGs_2rRNAs_NT_chain2 &
$phylobayes/tracecomp GTR_13PCGs_2rRNAs_NT_chain?.trace -o GTR_13PCGs_2rRNAs_NT_tracecomp
$phylobayes/bpcomp -x 764 10 GTR_13PCGs_2rRNAs_NT_chain1 GTR_13PCGs_2rRNAs_NT_chain2 -o GTR_13PCGs_2rRNAs_NT_bpcomp

# Dataset 3: 13PCGs_AA
mpirun -np 8 $phylobayes/pb_mpi -ncat 1 -gtr -d Data_21_Galeomorphii_13PCGs_AA.phy  GTR_AA_chain1 &
mpirun -np 8 $phylobayes/pb_mpi -ncat 1 -gtr -d Data_21_Galeomorphii_13PCGs_AA.phy GTR_AA_chain2 &
$phylobayes/tracecomp GTR_AA_chain?.trace -o GTR_AA_tracecomp
$phylobayes/bpcomp -x 1373 10 GTR_AA_chain1 GTR_AA_chain2 -o GTR_AA_bpcomp
```
2. The CAT-Poisson (CAT-F81) model
```
# Dataset 1: 13PCGs_NT
mpirun -np 24 $phylobayes/pb_mpi -cat -poisson -d Data_15_Galeomorphii_13PCGs_NT.phy CAT-POISSON_13PCGs_NT_chain1 &
mpirun -np 24 $phylobayes/pb_mpi -cat -poisson -d Data_15_Galeomorphii_13PCGs_NT.phy CAT-POISSON_13PCGs_NT_chain2 &
$phylobayes/tracecomp CAT-POISSON_13PCGs_NT_chain?.trace -o CAT-POISSON_13PCGs_NT_tracecomp
$phylobayes/bpcomp -x 2350 10 CAT-POISSON_13PCGs_NT_chain1 CAT-POISSON_13PCGs_NT_chain2 -o CAT-POISSON_13PCGs_NT_bpcomp

# Dataset 2: 13PCGs_2rRNAs_NT
mpirun -np 24 $phylobayes/pb_mpi -cat -poisson -d Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy CAT-POISSON_13PCGs_2rRNAs_NT_chain1 &
mpirun -np 24 $phylobayes/pb_mpi -cat -poisson -d Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy CAT-POISSON_13PCGs_2rRNAs_NT_chain2 &
$phylobayes/tracecomp CAT-POISSON_13PCGs_2rRNAs_NT_chain?.trace -o CAT-POISSON_13PCGs_2rRNAs_NT_tracecomp -x 3500
$phylobayes/bpcomp -x 3500 10 CAT-POISSON_13PCGs_2rRNAs_NT_chain1 CAT-POISSON_13PCGs_2rRNAs_NT_chain2 -o CAT-POISSON_13PCGs_2rRNAs_NT_bpcomp

# Dataset 3: 13PCGs_AA
mpirun -np 24 $phylobayes/pb_mpi -cat -poisson -d Data_21_Galeomorphii_13PCGs_AA.phy CAT-POISSON_AA_chain1 &
mpirun -np 24 $phylobayes/pb_mpi -cat -poisson -d Data_21_Galeomorphii_13PCGs_AA.phy CAT-POISSON_AA_chain2 &
$phylobayes/tracecomp CAT-POISSON_AA_chain?.trace -o CAT-POISSON_AA_tracecomp
$phylobayes/bpcomp -x 2601 10 CAT-POISSON_AA_chain1 CAT-POISSON_AA_chain2 -o CAT-POISSON_AA_bpcomp
```
3. The CAT-GTR model
```
# Dataset 1: 13PCGs_NT
mpirun -np 24 $phylobayes/pb_mpi -d Data_15_Galeomorphii_13PCGs_NT.phy CAT-GTR_13PCGs_NT_chain1
mpirun -np 24 $phylobayes/pb_mpi -d Data_15_Galeomorphii_13PCGs_NT.phy CAT-GTR_13PCGs_NT_chain2
$phylobayes/tracecomp CAT-GTR_13PCGs_NT_chain1.trace CAT-GTR_13PCGs_NT_chain2.trace -o CAT-GTR_13PCGs_NT_tracecomp
$phylobayes/bpcomp -x 7564 10 CAT-GTR_13PCGs_NT_chain1 CAT-GTR_13PCGs_NT_chain2 -o CAT-GTR_13PCGs_NT_bpcomp

# Dataset 2: 13PCGs_2rRNAs_NT
mpirun -np 24 $phylobayes/pb_mpi -d Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy CAT-GTR_13PCGs_2rRNAs_NT_chain1
mpirun -np 24 $phylobayes/pb_mpi -d Data_18_Galeomorphii_13PCGs_2rRNAs_NT.phy CAT-GTR_13PCGs_2rRNAs_NT_chain2
$phylobayes/tracecomp CAT-GTR_13PCGs_2rRNAs_NT_chain1.trace CAT-GTR_13PCGs_2rRNAs_NT_chain2.trace -o CAT-GTR_13PCGs_2rRNAs_NT_tracecomp
$phylobayes/bpcomp -x 6170 10 CAT-GTR_13PCGs_2rRNAs_NT_chain1 CAT-GTR_13PCGs_2rRNAs_NT_chain2 -o CAT-GTR_13PCGs_2rRNAs_NT_bpcomp

# Dataset 3: 13PCGs_AA
mpirun -np 24 $phylobayes/pb_mpi -d Data_21_Galeomorphii_13PCGs_AA.phy CAT-GTR_AA_chain1
mpirun -np 24 $phylobayes/pb_mpi -d Data_21_Galeomorphii_13PCGs_AA.phy CAT-GTR_AA_chain2
$phylobayes/tracecomp CAT-GTR_AA_chain1.trace CAT-GTR_AA_chain2.trace -o CAT-GTR_AA_tracecomp
$phylobayes/bpcomp -x 3900 10 CAT-GTR_AA_chain1 CAT-GTR_AA_chain2 -o CAT-GTR_AA_bpcomp
```
4. C60-GTR model (AA only)
```
# Dataset 3: 13PCGs_AA
mpirun -np 24 $phylobayes/pb_mpi -catfix c60 -gtr -d Data_21_Galeomorphii_13PCGs_AA.phy C60-GTR_AA_chain1
mpirun -np 24 $phylobayes/pb_mpi -catfix c60 -gtr -d Data_21_Galeomorphii_13PCGs_AA.phy C60-GTR_AA_chain2
$phylobayes/tracecomp C60-GTR_AA_chain?.trace -o C60-GTR_AA_tracecomp
$phylobayes/bpcomp -x 2469 10 C60-GTR_AA_chain1 C60-GTR_AA_chain2 -o C60-GTR_AA_bpcomp
```
### STEP 8.4: Check chain convergence

1. Generate site-specific log-likelihood statistics for the different models with the readpb_mpi program. Use a burn-in of 20% of the length of the MCMC chain and thin accordingly to obtain a total of at least 100 MCMC points given the computational intensiveness of estimating site log-likelihoods under the CAT models.
```
phylobayes = ./data

# Dataset 1: 13PCGs_NT
mpirun -np 48 $phylobayes/readpb_mpi -x 813 32 -sitelogl GTR_13PCGs_NT_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 813 32 -sitelogl GTR_13PCGs_NT_chain2 
mpirun -np 48 $phylobayes/readpb_mpi -x 2350 94 -sitelogl CAT-POISSON_13PCGs_NT_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 2350 94 -sitelogl CAT-POISSON_13PCGs_NT_chain2 
mpirun -np 48 $phylobayes/readpb_mpi -x 7564 300 -sitelogl CAT-GTR_13PCGs_NT_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 7564 300 -sitelogl CAT-GTR_13PCGs_NT_chain2

# Dataset 2: 13PCGs_rRNAs_NT
mpirun -np 48 $phylobayes/readpb_mpi -x 764 30 -sitelogl GTR_13PCGs_2rRNAs_NT_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 764 30 -sitelogl GTR_13PCGs_2rRNAs_NT_chain2 
mpirun -np 48 $phylobayes/readpb_mpi -x 3500 47 -sitelogl CAT-POISSON_13PCGs_2rRNAs_NT_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 3500 47 -sitelogl CAT-POISSON_13PCGs_2rRNAs_NT_chain2
mpirun -np 48 $phylobayes/readpb_mpi -x 6170 242 -sitelogl CAT-GTR_13PCGs_2rRNAs_NT_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 6170 242 -sitelogl CAT-GTR_13PCGs_2rRNAs_NT_chain2

# Dataset 3: 13PCGs_AA
mpirun -np 48 $phylobayes/readpb_mpi -x 1373 54 -sitelogl GTR_AA_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 1373 54 -sitelogl GTR_AA_chain
mpirun -np 48 $phylobayes/readpb_mpi -x 2469 98 -sitelogl C60-GTR_AA_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 2469 98 -sitelogl C60-GTR_AA_chain2
mpirun -np 48 $phylobayes/readpb_mpi -x 2601 104 -sitelogl CAT-POISSON_AA_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 2601 104 -sitelogl CAT-POISSON_AA_chain2
mpirun -np 48 $phylobayes/readpb_mpi -x 3900 154 -sitelogl CAT-GTR_AA_chain1
mpirun -np 48 $phylobayes/readpb_mpi -x 3900 154 -sitelogl CAT-GTR_AA_chain2
```
2. Collect the site-specific scores.
```
# Dataset 1: 13PCGs_NT
python3 ../scripts/read_loocv_waic.py GTR_13PCGs_NT_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py CAT-POISSON_13PCGs_NT_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py CAT-GTR_13PCGs_NT_chain?.sitelogl

# Dataset 2: 13PCGs_rRNAs_NT
python3 ../scripts/read_loocv_waic.py GTR_13PCGs_2rRNAs_NT_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py CAT-POISSON_13PCGs_2rRNAs_NT_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py CAT-GTR_13PCGs_2rRNAs_NT_chain?.sitelogl

# Dataset 3: 13PCGs_AA
python3 ../scripts/read_loocv_waic.py GTR_AA_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py CAT-POISSON_AA_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py C60-GTR_AA_chain?.sitelogl
python3 ../scripts/read_loocv_waic.py CAT-GTR_AA_chain?.sitelogl
```
3. Generate scores for the widely applicable information criterion (wAIC) and the asymptotically equivalent leave one-out cross-validation (LOO-CV) [7375] using the script read_loocv_waic.py, which is provided as part of the PHYLOBAYES_MPI package.
```
# Load the package
library(coda)

# Read the data
chain1 <- read.table("CAT-GTR_AA_chain1_reduced2matchchanin2.trace",header=TRUE)
chain2 <- read.table("CAT-GTR_AA_chain2.trace",header=TRUE)
chain3 <- read.table("CAT-GTR_13PCGs_NT_chain1.trace",header=TRUE)
chain4 <- read.table("CAT-GTR_13PCGs_NT_chain2_reduced2matchchanin1.trace",header=TRUE)
chain5 <- read.table("CAT-GTR_13PCGs_2rRNAs_NT_chain1.trace",header=TRUE)
chain6 <- read.table("CAT-GTR_13PCGs_2rRNAs_NT_chain2_reduced2matchchanin1.trace",header=TRUE)

# Create an MCMC object with the correct thinning interval (used 1000).
chain1 <- mcmc(chain1,thin=100)
chain2 <- mcmc(chain2,thin=100)
chain3 <- mcmc(chain3,thin=100)
chain4 <- mcmc(chain4,thin=100)
chain5 <- mcmc(chain5,thin=100)
chain6 <- mcmc(chain6,thin=100)

# Always start with a plot. A rough way of verifying convergence is to plot the trace and the posterior distribution of some of the parameters.
pdf("CAT-GTR_AA_chain1.trace.plot.pdf", width=10, height=10)
plot(chain1)
.=dev.off()
pdf("CAT-GTR_AA_chain2.trace.plot.pdf", width=10, height=10)
plot(chain2)
.=dev.off()
pdf("CAT-GTR_13PCGs_NT_chain1.trace.plot.pdf", width=10, height=10)
plot(chain3)
.=dev.off()
pdf("CAT-GTR_13PCGs_NT_chain2.trace.plot.pdf", width=10, height=10)
plot(chain5)
.=dev.off()
pdf("CAT-GTR_13PCGs_2rRNAs_NT_chain1.trace.plot.pdf", width=10, height=10)
plot(chain5)
.=dev.off()
pdf("CAT-GTR_13PCGs_2rRNAs_NT_chain2.trace.plot.pdf", width=10, height=10)
plot(chain6)
.=dev.off()

# Obtaining summary statistics for the MCMC chain. 
summary(chain1)
summary(chain2)
summary(chain3)
summary(chain4)
summary(chain5)
summary(chain6)

# Is the sample from the MCMC large enough? 
autocorr.diag(chain1)
effectiveSize(chain1)
autocorr.diag(chain2)
effectiveSize(chain2)
autocorr.diag(chain3)
effectiveSize(chain3)
autocorr.diag(chain4)
effectiveSize(chain4)
autocorr.diag(chain5)
effectiveSize(chain5)
autocorr.diag(chain6)
effectiveSize(chain6)

# Has the chain converged? 
combined.d1 = mcmc.list(chain1,chain2)
pdf("CAT-GTR_AA_chains-combined.trace.plot.pdf", width=10, height=10)
plot(combined.d1)
.=dev.off()
gelman.diag(combined.d1)
pdf("CAT-GTR_AA_chains-combined.trace.gelmanplot.pdf", width=10, height=10)
gelman.plot(combined.d1)
.=dev.off()
combined.d2 = mcmc.list(chain3,chain4)
pdf("CAT-GTR_13PCGs_NT_chains-combined.trace.plot.pdf", width=10, height=10)
plot(combined.d2)
.=dev.off()
gelman.diag(combined.d2)
pdf("CAT-GTR_13PCGs_NT_chains-combined.trace.gelmanplot.pdf", width=10, height=10)
gelman.plot(combined.d2)
.=dev.off()
combined.d3 = mcmc.list(chain5,chain6)
pdf("CAT-GTR_13PCGs_2rRNAs_NT_chains-combined.trace.plot.pdf", width=10, height=10)
plot(combined.d3)
.=dev.off()
gelman.diag(combined.d3)
pdf("CAT-GTR_13PCGs_2rRNAs_NT_chains-combined.trace.gelmanplot.pdf", width=10, height=10)
gelman.plot(combined.d3)
.=dev.off()
```
