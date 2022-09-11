# Bioinformatics pipeline of Winn et al. (2022) - Triakid Mitophlyogenomics
***[Simo: I am struggling to follow your workflow. May you please name your scripts appropriately i.e STEP-1_XXXX, STEP-2_XXX?]***

## STEP 1: Quality Control of Ion GeneStudio™ S5 data

***[Simo: How was quality control carried out? Provide a brief description of how the (1) sequence quality was checked, (2) low quality bases were trimmed and (3) adapters contamination was checked and removed.]***

## STEP 2: Mitogenome assemblies  

***[Simo: Please add a guideline for the graphical front end software and code for command line interphase software for the different assembly approaches you used]***

## STEP 3: Mitogenome annotations

***[Simo: Please add a guideline for the graphical front end software and code for command line interphase software for the different annotation approaches you used]***

## STEP 4: Sequence alignments and concatenation

### STEP 4.1: Ingroup and outgroup mitogenome retrieval from GenBank 
First, compile a list of mitochondrial genomes of the order Carcharhiniformes from Kousteni *et al*. (2021) and Wang *et al*. (2022), which include additional outgroup representatives from the orders Lamniformes and Orectolobiformes, and save the the text file as `kousteni-wang_mitogenomes_genbank.list`. Second, use the file `STEP-1_mitogenomes_genebak.list` to retrieve records from GenBank using [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) to download the entire mitogenome records in GenBank (full) format. Last, save the file as 'kousteni-wang.gb' in a folder named `1_data´, already containing newly assembled mitogenomes also in GenBank (full) format.

### STEP 4.2: Gene region extrations
First, within the `1_data´ folder, merge all the newely assembled mitogenome .gb files with those obtained from GenBank in the the preceeding step.
```
cat *.gb > winn_2022.gb
```
Second, extract rNRA and CDS DNA sequences from Genbank file using [GBSEQEXTRACTOR v.0.0.4](https://github.com/linzhi2013/gbseqextractor).
```
gbseqextractor -f winn_2022.gb -prefix winn_2022 -types rRNA -s # output file ´winn_2022.rrna.fasta´
gbseqextractor -f winn_2022.gb -prefix winn_2022 -types CDS -s	# output file `winn_2022.cds.fasta´
```
Third, merge the rNRA and CDS fasta files.
```
cat winn_2022.rrna.fasta winn_2022.cds.fasta > winn_2022.cds-rrna.fasta
	# Using Notepad+++, edit the file ´winn_2022.cds-rrna.fasta´ to standardize gene names.
	# For instance, some GenBank records call 12S rRNA as s-rRNA, 12S ribosomal RNA or rrnS,
	# CO1 as COX1, ND2 as nad2 etc. Standardize all genes to the following code: 
	# 'ATP6' 'ATP8' 'COX1' 'COX2' 'COX3' 'CYTB' 'ND1' 'ND2' 'ND3' 'ND4' 'ND4L' 'ND5' 'ND6' '12SrRNA' '16SrRNA'.
	# Save the edited file as `winn_2022.cds-rrna.std.fasta´.
```
Finally, extract individual gene sequences (.fa) from `winn_2022.cds-rrna.std.fasta´ using the custom script `maduna2022-gene-extractions.sh´:
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
Align protein coding genes using MACSE2 and rRNA genes using MAFFT (Q-INS-i algorithm).
First, we must create a new folder called `2_MSA´ and move .fa files of 13 PCGs to a subfolder called `13PCGs´ and the 2 rRNA genes to a subfolder called  `2RNAs´.
```
mkdir 13PCGs 2RNAs
mv *rRNA*.fa 2RNAs/
ls 2RNAs/*.fa | wc -l # folder should have 2 .fa files
mv *.fa 13PCGs/
ls 13PCGs/*.fa | wc -l # folder should have 13 .fa files
```
Second, download the the most recent release of [MACSE](https://bioweb.supagro.inra.fr/macse/) and copy the .jar file, 'macse_v2.06.jar' at present, into the ´13PCGs/´. We checked that sequecens were in the correct orientation and made adjustments accordingly in MEGA. Then we aligned the protein coding genes using the for loop script `maduna2022-13pcgs-msa.sh`.
```
datadir=./13PCGs
for i in $datadir/*.fa
do		
	java -jar -Xmx600m macse_v2.06.jar -prog alignSequences -seq "$i" -gc_def 2		
done
	# alignSequences: aligns nucleotide (NT) coding sequences using their amino acid (AA) translations.
	# gc_def: specify the genetic code 2 (The_Vertebrate_Mitochondrial_Code) or change accordingly for your study taxa.
```
Third, we opted to alig the rRNA genes using the online version of [MAFFT (version 7)](https://mafft.cbrc.jp/alignment/server/) using the Q-INS-i iterative refinement method, adjusting the direction according to the first sequence confirmed to be in the 5' -> 3' direction.
Fourth, we visually inspect and manually edit alignments using MEGA X. We removed terminal stop codons from the PCG alignments and trimmed alignements to a sequence lenght is divisible by 3. When needed, used BMGE (Block Mapping and Gathering with Entropy) to remove any remaining ambiguously aligned sites. For the RNAs, we used BMGE to remove ambiguously aligned sites.
Finally, before mitophylogenomic analysis, we produced three concatenated mitogenomic datasets from (i) the aligned individual PCGs datasets (Dataset 1: 13PCGs_NT dataset), (ii) the 13 PCGs plus the two rRNA genes (Dataset 2: 13PCGs_rRNAs_NT dataset) with the R package concatipede v1.0.1 (Vecchi and Bruneaux, 2021), and (iii) we derived the third mitogenomic dataset by translating the 13PCGs_NT dataset in MEGA (Dataset 3: 13PCGs_AA dataset). 

## STEP 5: Substitution saturation and data partitioning schemes 


## STEP 6: Phylogenetic reconstructions 



