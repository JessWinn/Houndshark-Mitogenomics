### Nucleotide Composition Plots ###

## Path to output PDFs
pdfPath <- 'C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/'

## Import the data

library(readxl)

C_amboinensis <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Carcharhinus amboinensis")
G_galeus <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Galeorhinus galeus")
G_melastomus <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Galeus melastomus")
H_japanica <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Hemitriakis japanica")
M_asterias <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Mustelus asterias")
M_griseus <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Mustelus griseus")
M_manazo <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Mustelus manazo")
M_mosis <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Mustelus mosis")
M_mustelus <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Mustelus mustelus")
M_palumbes <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Mustelus palumbes")
T_megalopterus <- read_excel("C:/Users/Jessica Winn/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/1_Mitogenome_Assembly/3_Genome_Annotation/2_Nucleotide Composition/Nucleotide Composition Plot/Nucleotide Composition.xlsx", sheet = "Triakis megalopterus")

## Initialize PDF
pdf(paste(pdfPath, '/Nucleotide_Composition.pdf', sep=''), width=8, height=6)

library(ggplot2)

## Carcharhinus ambionensis ##

ggplot(C_amboinensis) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title = "Carcharhinus amboinensis", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(C_amboinensis) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Carcharhinus amboinensis", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Galeorhinus galeus ##

ggplot(G_galeus) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Galeorhinus galeus", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(G_galeus) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Galeorhinus galeus", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Galeus melastomus ##

ggplot(G_melastomus) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Galeus melastomus", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(G_melastomus) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Galeus melastomus", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Hemitriakis japanica ##

ggplot(H_japanica) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Hemitriakis japanica", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(H_japanica) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Hemitriakis japanica", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Mustelus asterias ##

ggplot(M_asterias) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus asterias", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(M_asterias) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus asterias", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Mustelus griseus ##

ggplot(M_griseus) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus griseus", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(M_griseus) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus griseus", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Mustelus manazo ##

ggplot(M_manazo) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus manazo", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(M_manazo) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus manazo", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Mustelus mosis ##

ggplot(M_mosis) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus mosis", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(M_mosis) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus mosis", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Mustelus mustelus ##

ggplot(M_mustelus) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus mustelus", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(M_mustelus) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus mustelus", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Mustelus palumbes ##

ggplot(M_palumbes) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus palumbes", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(M_palumbes) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Mustelus palumbes", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

## Triakis megalopterus ##

ggplot(T_megalopterus) + 
  geom_line(aes(x=Gene, y=AT_content, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_content, group=1), color = "darkturquoise", size=1) +
  labs(title="Triakis megalopterus", x="Gene", y="Content(%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

ggplot(T_megalopterus) + 
  geom_line(aes(x=Gene, y=AT_skew, group=1), color = "coral", size=1) + 
  geom_line(aes(x=Gene, y=GC_skew, group=1), color = "darkturquoise", size=1) +
  labs(title="Triakis megalopterus", x="Gene", y="Skewness") +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=18, face='italic')) +
  theme(axis.text.x = element_text(angle = 52, size=9))

.=dev.off()