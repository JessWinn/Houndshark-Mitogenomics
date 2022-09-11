# Bioinformatics pipeline of Winn et al. (2022) - Triakid Mitophlyogenomics
***[Simo: I am struggling to follow your workflow. May you please name your scripts appropriately i.e STEP-1_XXXX, STEP-2_XXX?]***

## STEP 1: Quality Control of Ion GeneStudio™ S5 data
***[Simo: How was quality control carried out? Provide a brief description of how the (1) sequence quality was checked, (2) low quality bases were trimmed and (3) adapters contamination was checked and removed.]***

## STEP 2: Mitogenome assemblies  

## STEP 3: Mitogenome annotations

## STEP 4: Mitophylogenomics

### STEP 4.1: Ingroup and outgroup mitogenome retrieval from GenBank 
First, compile a list of mitochondrial genomes of the order Carcharhiniformes from Kousteni *et al*. (2021) and Wang *et al*. (2022), which include additional outgroup representatives from the orders Lamniformes and Orectolobiformes, and save the the text file as ***STEP-1_mitogenomes_genebak.list***. Second, use the file ***STEP-1_mitogenomes_genebak.list*** to retrieve records from GenBank using [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) to download the entire mitogenome records in GenBank (full) format. Last, save the file as 'kousteni-wang.gb' in a folder named ***data***, already containing newly assembled mitogenomes also in GenBank (full) format.

### STEP 4.2: Gene region extrations
First, within the ***data*** folder, merge all the newely assembled mitogenome .gb files with those obtained from GenBank in the the preceeding step.
```
cat *.gb > winn_2022.gb
```
Second, extract rNRA and CDS DNA sequences from Genbank file using [GBSEQEXTRACTOR v.0.0.4](https://github.com/linzhi2013/gbseqextractor)
```
gbseqextractor -f winn_2022.gb -prefix winn_2022 -types rRNA -s # output file 'winn_2022.rrna.fasta'
gbseqextractor -f winn_2022.gb -prefix winn_2022 -types CDS -s	# output file 'winn_2022.cds.fasta'
```


