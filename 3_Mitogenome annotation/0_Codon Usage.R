### Codon Usage Plots ###

## Carcharhinus ambionensis ##

## Path to output PDFs
pdfPath <- 'C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform mitogenomics/Codon Usage'

# Import the data

library(readxl)

C_ambionensis <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Carcharhinus ambionensis")
G_galeus <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Galeorhinus galeus")
G_melastomus <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Galeus melastomus")
H_japanica <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Hemitriakis japanica")
M_asterias <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Mustelus asterias")
M_griseus <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Mustelus griseus")
M_manazo <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Mustelus manazo")
M_mosis <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Mustelus mosis")
M_mustelus <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Mustelus mustelus")
M_palumbes <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Mustelus palumbes")
T_megalopterus <- read_excel("C:/Users/peter.winn/OneDrive/Jess OneDrive/OneDrive/Documents/Masters/Ion Torrent/Carcharhiniform Mitogenomics/Codon Usage/Codon Usage.xlsx", sheet = "Triakis megalopterus")

# Amino Acid composition

### Initialize PDF
pdf(paste(pdfPath, '/AA_composition.pdf', sep=''), width=8, height=6)

library(ggplot2)

#C_ambionensis

ggplot(C_ambionensis, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(C_ambionensis, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#G_galeus

ggplot(G_galeus, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(G_galeus, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#G_melastomus

ggplot(G_melastomus, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(G_melastomus, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#H_japanica

ggplot(H_japanica, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(H_japanica, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#M_asterias

ggplot(M_asterias, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(M_asterias, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#M_griseus

ggplot(M_griseus, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(M_griseus, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#M_manazo

ggplot(M_manazo, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(M_manazo, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#M_mosis

ggplot(M_mosis, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(M_mosis, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#M_mustelus

ggplot(M_mustelus, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(M_mustelus, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#M_palumbes

ggplot(M_palumbes, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(M_palumbes, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

#T_megalopterus

ggplot(T_megalopterus, aes(x = Amino_Acid, y = N)) +
  geom_bar(stat = "identity", fill = 4) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic'))

ggplot(T_megalopterus, aes(x = Amino_Acid, y=RSCU, fill=Codon, label=Codon2)) +
  geom_bar(stat = "identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5, size=20, face='italic')) +
  theme(legend.position="none") +
  scale_fill_manual('Codon', values=c('#FFA54F', '#00C1AA', '#FF69B4', '#9370DB'))

.=dev.off()