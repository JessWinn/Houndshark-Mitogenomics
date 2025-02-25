### STEP 6.5

## Adapted from Lanfear’s R script (Minh et al., 2020).

library(viridis)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(GGally)
library(entropy)

pdfPath <- './CF plots/'
pdf(paste(pdfPath, '/CF_plots.pdf', sep=''), width=8, height=6)

# PS1

PS1 = read.delim("./PS1_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS1)[18] = "bootstrap"
names(PS1)[19] = "branchlength"


# plot the values
ggplot(PS1, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS1AA

PS1AA = read.delim("./PS1AA_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS1AA)[18] = "bootstrap"
names(PS1AA)[19] = "branchlength"


# plot the values
ggplot(PS1AA, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS2

PS2 = read.delim("./PS2_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS2)[18] = "bootstrap"
names(PS2)[19] = "branchlength"


# plot the values
ggplot(PS2, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS3

PS3 = read.delim("./PS3_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS3)[18] = "bootstrap"
names(PS3)[19] = "branchlength"


# plot the values
ggplot(PS3, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS4

PS4 = read.delim("./PS4_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS4)[18] = "bootstrap"
names(PS4)[19] = "branchlength"


# plot the values
ggplot(PS4, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS5

PS5 = read.delim("./PS5_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS5)[18] = "bootstrap"
names(PS5)[19] = "branchlength"


# plot the values
ggplot(PS5, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS5AA

PS5AA = read.delim("./PS5AA_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS5AA)[18] = "bootstrap"
names(PS5AA)[19] = "branchlength"


# plot the values
ggplot(PS5AA, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS6

PS6 = read.delim("./PS6_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS6)[18] = "bootstrap"
names(PS6)[19] = "branchlength"


# plot the values
ggplot(PS6, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS7

PS7 = read.delim("./PS7_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS7)[18] = "bootstrap"
names(PS7)[19] = "branchlength"


# plot the values
ggplot(PS7, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# PS8

PS8 = read.delim("./PS8_run03_concord.cf.stat", header = T, comment.char = '#')

names(PS8)[18] = "bootstrap"
names(PS8)[19] = "branchlength"


# plot the values
ggplot(PS8, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) + 
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

.=dev.off()