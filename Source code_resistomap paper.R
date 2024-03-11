#Fig 1 Number of detected genes
getwd()

library(readr)
f2numgene <- read_csv("/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon/f2numgene.csv", 
                      col_types = cols(Aminoglycoside = col_number(), 
                                       `Beta Lactam` = col_number(), Vancomycin = col_number(), 
                                       Other = col_number(), MDR = col_number(), 
                                       Sulfonamide = col_number(), Tetracycline = col_number(), 
                                       Trimethoprim = col_number(), Phenicol = col_number(), 
                                       Integrons = col_number(), MGE = col_number(), 
                                       Quinolone = col_number(), MLSB = col_number()))
View(f2numgene)

head(f2numgene)
# Count the number of columns
f2numgenen <- ncol(f2numgene)
f2numgenel <- nrow(f2numgene)
# Print the number of columns
print(f2numgenel)

# Install and load the reshape2 package
library(reshape2)

# Assuming your original data frame is named 'original_df'

# Melt the data frame to convert it from wide format to long format
melted_f2numgene <- melt(f2numgene, id.vars = c("Sample number", "Description", "Site"))
print(melted_f2numgene)

# Create data frame
numgene <- as.data.frame(melted_f2numgene)
print(numgene)

# Install and load the openxlsx package
library(openxlsx)

# Assuming your data frame is named 'data_frame' and you want to save it to 'output_file.xlsx'

# Specify the file path for the output Excel file
output_file <- "numgene_file.xlsx"

# Save the data frame to the Excel file
write.xlsx(numgene, file = output_file, rowNames = FALSE)

# Print a message indicating the successful save
cat("Data frame saved to Excel file:", output_file, "\n")

# import data for stack bar
library(readr)
numgene <- read_csv("numgene.csv", col_types = cols(count = col_number()))
View(numgene)
head(numgene)

ID <- as.factor(numgene$ID)
count <- sample(numgene$count)
class <- sample(numgene$Class)

df <- data.frame(x = numgene$ID, y = numgene$count, group = numgene$Class)
head(df)
View(df)

# install.packages("ggplot2")
library(ggplot2)
library(tidyverse)

p1 <- df %>%
  ggplot(aes(x = x, y = y, fill = group, label = y)) +
  geom_bar(stat = "identity", position = "stack", width=0.6, colour="black", linewidth=0.2) +
  facet_grid(~factor(numgene$Site, levels=c('Wastewater', 'Aeration tank', 'Effluent', 'Downstream', 'Upstream')), scales = "free_x") + 
  theme(axis.title.x = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", family = "Arial", size = 7), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #site detail label
        legend.key.height= unit(0.2, 'cm'),
        legend.key.width= unit(0.2, 'cm'),
        axis.text = element_text(color = "black", family = "Arial", size = 5),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size= 0.5)) + xlab("Sampling sites") + 
  ylab("Number of detected genes") + labs(fill = "Antibiotics Class") +
  scale_fill_manual(values = colorRampPalette(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"))(13))
p1
ggsave(file="Stakedbar2.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 190, height = 80, units = "mm", dpi=600) #190 #240

# Fig 2AB Heat map Relatives abundance to 16S rRNA gene
library(ggplot2)
library(rstatix)
library(reshape2)
library(ggpubr)
library(dplyr)
library(survival)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)

# data with group (2A)
library(readxl)
heatmap_gr2_1 <- read_excel("heatmap_gr2.1.xlsx", 
                            col_types = c("text", "text", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric"))
View(heatmap_gr2_1)


head(heatmap_gr2_1)

# Count the number of columns
heatmap_gr2_1col <- ncol(heatmap_gr2_1)
# Print the number of columns #20
print(heatmap_gr2_1col)

####reshape red heatmap
# cut data frame cut group out
column1 <- as.data.frame(heatmap_gr2_1[,3:20])
print(column1)

#add columns to data frame
column1$gene <- heatmap_gr2_1$gene
head(column1)

column1melt <- melt(column1) 
head(column1melt) #see results only head
column1melt <- as.data.frame(column1melt)

ggp <- ggplot(column1melt, aes(variable, gene)) +
  geom_tile(aes(fill = value), colour = "white") +
  facet_grid(heatmap_gr2_1$group~., scales = "free_y", space = "free_y")

ggp + scale_fill_gradient2(low = "white", mid = "white", high = "blue4") +
  xlab("Sampling sites") + 
  ylab("ARGs and MGEs") + labs(fill = "Relative Abundance") +
  theme(axis.title.x = element_text(color = "black", size = 5, face = "bold"), #ชื่อแกน
        axis.title.y = element_text(color = "black", size = 5, face = "bold"), 
        legend.text = element_text(color = "black", size = 5), #Inlet Outlet #criteria label
        legend.title = element_text(color = "black", size = 5, face = "bold"), #Sampling point #criteria label
        strip.text.y = element_text(color = "black", size = 4, face = "bold"),
        axis.text.y = element_text(color = "black", size = 5, face = "italic"), #y label
        axis.text.x = element_text(color = "black", size = 5) #x label
  )

ggsave(file="heatmap_resistomegr1.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 150, height = 160, units = "mm", dpi=600)

# data with group (2B)
library(readxl)
heatmap_gr1_1 <- read_excel("heatmap_gr1.1.xlsx", 
                            col_types = c("text", "text", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric"))
View(heatmap_gr1_1)

# Count the number of columns
heatmap_gr1_1col <- ncol(heatmap_gr1_1)
# Print the number of columns #20
print(heatmap_gr1_1col)

####reshape red heatmap
# cut data frame cut group out
column2 <- as.data.frame(heatmap_gr1_1[,3:20])
print(column2)

#add columns to data frame
column2$gene <- heatmap_gr1_1$gene
head(column2)

column2melt <- melt(column2) 
head(column2melt) #see results only head
column2melt <- as.data.frame(column2melt)

ggp <- ggplot(column2melt, aes(variable, gene)) +
  geom_tile(aes(fill = value), colour = "white") +
  facet_grid(heatmap_gr1_1$group~., scales = "free_y", space = "free_y")

ggp + scale_fill_gradient2(low = "white", mid = "white", high = "blue4") +
  xlab("Sampling sites") + 
  ylab("ARGs and MGEs") + labs(fill = "Relative Abundance") +
  theme(axis.title.x = element_text(color = "black", size = 5, face = "bold"), #ชื่อแกน
        axis.title.y = element_text(color = "black", size = 5, face = "bold"), 
        legend.text = element_text(color = "black", size = 5), #Inlet Outlet #criteria label
        legend.title = element_text(color = "black", size = 5, face = "bold"), #Sampling point #criteria label
        strip.text.y = element_text(color = "black", size = 4, face = "bold"),
        axis.text.y = element_text(color = "black", size = 5, face = "italic"), #y label
        axis.text.x = element_text(color = "black", size = 5) #x label
  )

ggsave(file="heatmap_resistomegr2.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 150, height = 160, units = "mm", dpi=600)

## Figure 3 NDMS
library(ggplot2)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)
library(tidyverse)
library(permute)
library(lattice)
library(vegan)
library(ggvegan)
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)


#1. import data
library(readr)
NMDS3 <- read_csv("NMDS3.csv")
View(NMDS3)
# Count the number of columns
num_cols3 <- ncol(NMDS3)
# Print the number of columns #98
print(num_cols3)
# Convert all 98 columns to numeric
data3 <- NMDS3[,4:98]
print(data3)
data3 <- as.data.frame(lapply(data3, as.numeric), stringsAsFactors = FALSE)
print(data3)
# Check for NA values in each column #0
sapply(data3, function(x) sum(is.na(x)))

#2. dissimilarities matrix calculation (2 type, more common = bray, less common = helinger)
#bray methods
data3_dissmatrixbray <- vegdist(data3, method = "bray")
print(data3_dissmatrixbray)
dis.mat<-sqrt(data3_dissmatrixbray) #no run

#3. metaMDS calculation
#3.3 data 3 #0.08639729
data3_result4 <- metaMDS(data3_dissmatrixbray, distance='bray', species = species_data)
####plot graph stress  # Shepards test/goodness of fit
stressplot(data3_result4) 
goodness(data3_result4)
####stress value 0.05-0.1=good fit (#stress 0.05) #0.08639729
data3_result4$stress
###Print the NMDS result
print(data3_result4)
head(data3_result4)

#4. store diss matrix for plot graph
data.scores3 <- as.data.frame(scores(data3_result4))
print(data.scores3)
#add columns to data frame
data.scores3$assay <- NMDS3$assay
data.scores3$site <- NMDS3$site
data.scores3$season <- NMDS3$season
head(data.scores3)

#5. Plot output
jpeg(file="NDMSKa.jpeg", width= 190, height = 150, units = "mm", bg = "white", res = 600) #190 #240
A <- ordiplot(data3_result4, type = "n", xlim=c(-1, 0.7), ylim = c(-1, 0.5))

B <- nmds_ordiellipse <- cbind(data3_result4, data.scores3$site)

cols <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2')

C <- ordiellipse(data3_result4, data.scores3$site, col = cols,
                 label = FALSE, kind = "se", conf=0.95,
                 lwd = 2, lty = "dashed")
D <- ordispider(data3_result4, groups = data.scores3$site,
                col = cols, kind = "se", label = TRUE, lwd = 2)
E <- points(data3_result4, groups = data.scores3$site, pch=24, bg = 'yellow') #same same

E1 <- legend('bottomright', legend = "Bray Curtis, stress = 0.086", 
             inset = c(0.02, 0.02), bty = "n")
E2 <- legend(x="bottomleft", legend=c("AE: aeration tank", "CL2: effluent after disinfection", 
                                      "DW: downstream canal", "UP: upstream canal", 
                                      "WW: untreated wastewater"),
             col=cols, lwd = 2, inset = c(0.02, 0.02), bty = "n")
dev.off()

#6. #Permanova
#input data to matrix
library(readr)
NMDS3 <- read_csv("NMDS3.csv")
View(NMDS3)

data3 <- NMDS3[,4:98]
data3.matrix <- as.matrix(data3)
NMDS3$site <- as.factor(NMDS3$site)
head(data3.matrix)

# NMDS3$site R2 = 4.1987, p-value = 0.001
data3.mat <- sqrt(data3.matrix)
data3.dist <- vegdist(data3.mat, method='bray')

# multilevel pairwise
data3.div <- adonis2(data3.dist~NMDS3$site, data=NMDS3, permutations = 999, method="bray")
data3.div

##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

## end copy here

data3.divp <- pairwise.adonis(NMDS3[,4:98], NMDS3$site)
data3.divp

# Figure 4 boxplot qpcr with sig
library(tidyverse)
library(gapminder)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(corrplot)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)
library(ggsignif)
library(readxl)
library(readxl)
library(patchwork)
library(cowplot)
library(lattice)
library(PASWR)

X1INT <- read_excel("1INT.xlsx", col_types = c("text", 
                                               "text", "numeric"))
X2AE <- read_excel("2AE.xlsx", col_types = c("text", 
                                             "text", "numeric"))
X3TEM <- read_excel("3TEM.xlsx", col_types = c("text", 
                                               "text", "numeric"))
X4NDM <- read_excel("4NDM.xlsx", col_types = c("text", 
                                               "text", "numeric"))
X5VIM <- read_excel("5VIM.xlsx", col_types = c("text", 
                                               "text", "numeric"))
X6tet <- read_excel("6tet.xlsx", col_types = c("text", 
                                               "text", "numeric"))
X7mcr <- read_excel("7mcr.xlsx", col_types = c("text", 
                                               "text", "numeric"))
X8cr <- read_excel("8cr.xlsx", col_types = c("text", 
                                             "text", "numeric"))
X9sRNA <- read_excel("9sRNA.xlsx", col_types = c("text", 
                                                 "text", "numeric"))



p1 <- ggplot(X1INT, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(10.1, 10.3, 10.5)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p1
save_plot("p1.jpeg", p1)

p2 <- ggplot(X2AE, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p2
save_plot("p2.jpeg", p2)

p3 <- ggplot(X3TEM, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(9.1, 9.4, 9.7)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p3
save_plot("p3.jpeg", p3)

p4 <- ggplot(X4NDM, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(9.1, 9.4, 9.7)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p4
save_plot("p3.jpeg", p4)

p5 <- ggplot(X5VIM, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("CL2", "UP"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(9.1, 9.4, 9.7, 10)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p5
save_plot("p5.jpeg", p5)


p6 <- ggplot(X6tet, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(8.1, 8.4, 8.7, 9)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p6
save_plot("p6.jpeg", p6)

p7 <- ggplot(X7mcr, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(4.9, 5.2, 5.5, 5.8)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p7
save_plot("p7.jpeg", p7)

p8 <- ggplot(X8cr, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "AE"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(6, 6.3)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("Sampling sites") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p8
save_plot("p8.jpeg", p8)

p9 <- ggplot(X9sRNA, aes(x=Site, y=TC, fill =Site)) +
  geom_boxplot(alpha=1, width=0.4, size=0.1) + facet_wrap(~Name, ncol = 2, scales = "free_y") +
  geom_signif(comparisons = list(c("WW", "CL2"), c("WW", "DW"), c("WW", "UP")),
              map_signif_level = TRUE, textsize = 3,
              y_position = c(9.7, 10, 10.3)) + 
  stat_boxplot(geom= 'errorbar' , width = 0.4, alpha=1, position = position_dodge(width = 1)) +
  xlab("Sampling sites") + labs(fill = "Sampling site") + ylab("Log10 copies/100 mL") +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP")) +
  scale_fill_manual(values=c("#ABE188", "#F7EF99", "#F1BB87", "#F78E69", "#5D675B"))
p9
save_plot("p9.jpeg", p9)

fig1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 2,
                  labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"), label_size = 14)

ggsave(file="9figboxp.jpeg", fig1, path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 180, height = 235, units = "mm", dpi=600)

# Figure 5 boxplot removal 
library(readxl)
boxplot_remove <- read_excel("boxplot_remove.xlsx", 
                             col_types = c("text", "text", "numeric"))
View(boxplot_remove)

boxplot_remove %>% 
  as_tibble() %>% 
  ggplot(aes(x=Site, y=MARI, fill =Site)) + 
  geom_boxplot(alpha=1, width=0.4, size=0.1) + 
  stat_boxplot(geom= 'errorbar', width = 0.2, alpha=1, position = position_dodge(width = 0.4)) +
  xlab("ARGs and MGEs") + 
  ylab("log10 reduction values (LRV)") + labs(fill = "ARGs and MGEs") +
  theme(axis.title.x = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", family = "Arial", size = 7), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #site detail label
        plot.margin = margin(1, 1, 1, 1, "mm"),
        axis.text.x = element_text(color = "black", family = "Arial", size = 7),
        axis.text.y = element_text(color = "black", family = "Arial", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  scale_fill_manual(values = colorRampPalette(c("#F9ED69", "#F08A5D", "#B83B5E", "#6A2C70"))(8))

ggsave(file="boxplotremove.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 90, height = 90, units = "mm", dpi=600)

# Figure 6 qPCR WQ #absolute
library(ggcorrplot)

library(readr)
qPCRab_WQ <- read_csv("qPCRab_WQ.csv", col_types = cols(DO = col_number(), 
                                                        BOD = col_number(), COD = col_number(), 
                                                        TSS = col_number(), TDS = col_number(), 
                                                        Temp = col_number(), pH = col_number(), 
                                                        Salinity = col_number(), Turbidity = col_number(), 
                                                        INT = col_number(), SUL = col_number(), 
                                                        TEM = col_number(), NDM = col_number(), 
                                                        VIM = col_number(), TETQ = col_number(), 
                                                        MCR = col_number(), CrAss = col_number()))
View(qPCRab_WQ)
attach(qPCRab_WQ)
head(qPCRab_WQ)
cor(qPCRab_WQ)
p.mat <- cor_pmat(qPCRab_WQ)
head(p.mat)
m <- cor(qPCRab_WQ)

p<- ggcorrplot(m, hc.order = F, type = "lower",
               p.mat = p.mat, sig.level = .05, tl.cex = 7, pch = NA, insig = c("blank"), pch.cex = NA,
               lab = TRUE, outline.color = "white", lab_col = "black", lab_size = 2,
               ggtheme = ggplot2::theme_gray, colors = c("#FFCC70", "#FFFADD", "#8ECDDD")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 4), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 4, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 4),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.2)) 
p

ggsave(file="matrixqPCR_WQ3.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 100, height = 90, units = "mm", dpi=600)


# new corr plot
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)

corr <- round(cor(qPCRab_WQ), 1)

p.df <- as.data.frame(ggcorrplot::cor_pmat(qPCRab_WQ))

labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#FFCC70", "#FFFADD", "#8ECDDD")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 2),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

cor.plot.labs
ggsave(file="matrixqPCR_WQ4.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 200, height = 190, units = "mm", dpi=600)


# Figure 7 correlation phage and qPCR arg only up and down
# without MCR because the standard deviation is zero
library(readr)
cor_phageARG <- read_csv("cor_phageARG.csv", 
                         col_types = cols(CrAss = col_number(), 
                                          INT = col_number(),NDM = col_number(), SUL = col_number(), 
                                          TEM = col_number(), TETQ = col_number(), 
                                          VIM = col_number()))
attach(cor_phageARG)
head(cor_phageARG)

p.mat <- cor_pmat(cor_phageARG) #correlation matrix with p-values
head(p.mat)
p.mat
m <- cor(cor_phageARG) 
m

ggcorrplot(m, hc.order = TRUE, type = "lower", 
           p.mat = p.mat, sig.level = .05, tl.cex = 7, pch = 4, insig = c("pch"), pch.cex = 3,
           lab = TRUE, outline.color = "white", lab_col = "black", lab_size = 3,
           ggtheme = ggplot2::theme_gray, colors = c("#F4EEEE", "#FFDBAA", "#96C291")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 4), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 4, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 4),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.2))

ggsave(file="matrixARG_phage.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 80, height = 70, units = "mm", dpi=600)

# new corr plot
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)

corr <- round(cor(cor_phageARG), 1)

p.df <- as.data.frame(ggcorrplot::cor_pmat(cor_phageARG))

labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#F4EEEE", "#FFDBAA", "#96C291")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 2),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

cor.plot.labs
ggsave(file="matrixqPCR_WQ56.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 100, height = 90, units = "mm", dpi=600)

#Figure 8 supplement Ht-qpcr and qPCR
library(ggplot2)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)

library(readxl)
HTQ_data <- read_excel("HTQ_data.xlsx", col_types = c("text", 
                                                      "text", "numeric", "numeric"))
View(HTQ_data)
head(HTQ_data)

#row vs name
HTQ_data_INT <- HTQ_data[1:18, ]
HTQ_data_SUL <- HTQ_data[19:36, ]
HTQ_data_TEM <- HTQ_data[37:54, ]
HTQ_data_NDM <- HTQ_data[55:72, ]
HTQ_data_VIM<- HTQ_data[73:90, ]
HTQ_data_TETQ <- HTQ_data[91:108, ]
HTQ_data_MCR <- HTQ_data[109:126, ]
HTQ_data_MCR <- HTQ_data[127:144, ]
head(HTQ_data_INT)
head(HTQ_data_MCR)

HTQ_data %>% 
  as_tibble() %>%
  gather(wet, Concentrations, -Name, -Site) %>% 
  ggplot(aes(Site, Concentrations, fill =wet)) + 
  geom_boxplot(alpha=1, width=0.5, size=0.1) + 
  facet_wrap(~Name, ncol = 2, scales = "free_y") + 
  stat_boxplot(geom= 'errorbar' , width = 0.25, alpha=1, position = position_dodge(width = 0.5)) +
  xlab("Sampling sites") + 
  ylab("Relative abundance (copies/16s rRNA)") + labs(fill = "Methods") + 
  theme(axis.title.x = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", family = "Arial", size = 7), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 7, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 7),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1)) + 
  scale_x_discrete(limits=c("WW", "AE", "CL2", "DW", "UP"))

ggsave(file="boxplotSite.jpeg", path = "/Users/ydpr/Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon",
       width= 185, height = 235, units = "mm", dpi=600) #190 #240