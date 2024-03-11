#Import Data

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

###################Methods and Season_sum qPCR and ht-PCR##########################

library(readr)
Qseason_all <- read_csv("Library/CloudStorage/OneDrive-WalailakUniversity/MESTH_OneDrive/ปี 1 เทอม 3/Thesis/Five Chapter/บท 4/Python and R/R/Thesis Analysis/Stat_NFon/Qseason_all.csv", 
                        col_types = cols(INT_q = col_number(), 
                                         SUL_q = col_number(), TEM_q = col_number(), 
                                         NDM_q = col_number(), VIM_q = col_number(), 
                                         TETQ_q = col_number(), MCR_q = col_number(), 
                                         crAssphage = col_number()))
View(Qseason_all)
head(Qseason_all)

Dry <- Qseason_all[1:12,]
Wet <- Qseason_all[13:24,]

Qseason_all$Season <- as.factor(Qseason_all$Season)

sum.stat <- cenfit(Qseason_all$MCR_q, Qseason_all$MCR_cen, Qseason_all$Season) 
sum.stat

########NDM###########
summary(Dry$NDM_q)
sd(Dry$NDM_q)
shapiro.test(Dry$NDM_q)
summary(Wet$NDM_q)
sd(Wet$NDM_q)
shapiro.test(Wet$NDM_q)

########TEM###########
summary(Dry$TEM_q)
sd(Dry$TEM_q)
shapiro.test(Dry$TEM_q)
summary(Wet$TEM_q)
sd(Wet$TEM_q)
shapiro.test(Wet$TEM_q)

########VIM###########
summary(Dry$VIM_q)
sd(Dry$VIM_q)
shapiro.test(Dry$VIM_q)
summary(Wet$VIM_q)
sd(Wet$VIM_q)
shapiro.test(Wet$VIM_q)

########INT###########
summary(Dry$INT_q)
sd(Dry$INT_q)
shapiro.test(Dry$INT_q)
summary(Wet$INT_q)
sd(Wet$INT_q)
shapiro.test(Wet$INT_q)

########MCR###########
summary(Dry$MCR_q)
sd(Dry$MCR_q)
shapiro.test(Dry$MCR_q)
summary(Wet$MCR_q)
sd(Wet$MCR_q)
shapiro.test(Wet$MCR_q)

########SUL###########
summary(Dry$SUL_q)
sd(Dry$SUL_q)
shapiro.test(Dry$SUL_q)
summary(Wet$SUL_q)
sd(Wet$SUL_q)
shapiro.test(Wet$SUL_q)

########TETQ###########
summary(Dry$TETQ_q)
sd(Dry$TETQ_q)
shapiro.test(Dry$TETQ_q)
summary(Wet$TETQ_q)
sd(Wet$TETQ_q)
shapiro.test(Wet$TETQ_q)


library(readr)
Hseason <- read_csv("Hseason.csv", col_types = cols(INT_h = col_number(), 
                                                    SUL_h = col_number(), TEM_h = col_number(), 
                                                    NDM_h = col_number(), VIM_h = col_number(), 
                                                    TETQ_h = col_number(), MCR_h = col_number()))


head(Hseason)

Dry <- Hseason[1:9,]
Wet <- Hseason[10:18,]

Hseason$Season <- as.factor(Hseason$Season)
sum.stat <- cenfit(Hseason$MCR_h, Hseason$MCR_cen, Hseason$Season) 
sum.stat

########Phage###########
summary(Dry$Phage)
sd(Dry$Phage)
shapiro.test(Dry$Phage)
summary(Wet$Phage)
sd(Wet$Phage)
shapiro.test(Wet$Phage)

########NDM###########
summary(Dry$NDM_h)
sd(Dry$NDM_h)
shapiro.test(Dry$NDM_h)
summary(Wet$NDM_h)
sd(Wet$NDM_h)
shapiro.test(Wet$NDM_h)

########TEM###########
summary(Dry$TEM_h)
sd(Dry$TEM_h)
shapiro.test(Dry$TEM_h)
summary(Wet$TEM_h)
sd(Wet$TEM_h)
shapiro.test(Wet$TEM_h)

########VIM###########
summary(Dry$VIM_h)
sd(Dry$VIM_h)
shapiro.test(Dry$VIM_h)
summary(Wet$VIM_h)
sd(Wet$VIM_h)
shapiro.test(Wet$VIM_h)

########INT###########
summary(Dry$INT_h)
sd(Dry$INT_h)
shapiro.test(Dry$INT_h)
summary(Wet$INT_h)
sd(Wet$INT_h)
shapiro.test(Wet$INT_h)

########MCR###########
summary(Dry$MCR_h)
sd(Dry$MCR_h)
shapiro.test(Dry$MCR_h)
summary(Wet$MCR_h)
sd(Wet$MCR_h)
shapiro.test(Wet$MCR_h)

########SUL###########
summary(Dry$SUL_h)
sd(Dry$SUL_h)
shapiro.test(Dry$SUL_h)
summary(Wet$SUL_h)
sd(Wet$SUL_h)
shapiro.test(Wet$SUL_h)

########TETQ###########
summary(Dry$TETQ_h)
sd(Dry$TETQ_h)
shapiro.test(Dry$TETQ_h)
summary(Wet$TETQ_h)
sd(Wet$TETQ_h)
shapiro.test(Wet$TETQ_h)


#####Season diff########
library(readr)
Qseason <- read_csv("Qseason.csv", col_types = cols(INT_q = col_number(), 
                                                    SUL_q = col_number(), TEM_q = col_number(), 
                                                    NDM_q = col_number(), VIM_q = col_number(), 
                                                    TETQ_q = col_number(), MCR_q = col_number()))

head(Qseason)
Dry <- Qseason[1:9,]
Wet <- Qseason[10:18,]

###########NDM#############
shapiro.test(Dry$NDM_q) #0.06522
shapiro.test(Wet$NDM_q) #0.02035
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$NDM_q, Wet$NDM_q)                                        #0.6665
wilcox.test(Dry$NDM_q, Wet$NDM_q, alternative="greater")                 #0.3332
wilcox.test(Dry$NDM_q, Wet$NDM_q, alternative="less")                    #0.6976
###########TEM######################
shapiro.test(Dry$TEM_q) #0.05853
shapiro.test(Wet$TEM_q) #0.01296
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$TEM_q, Wet$TEM_q)                                        
wilcox.test(Dry$TEM_q, Wet$TEM_q, alternative="greater")
wilcox.test(Dry$TEM_q, Wet$TEM_q, alternative="less")  
###########VIM######################
shapiro.test(Dry$VIM_q) #0.2347
shapiro.test(Wet$VIM_q) #0.004783
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$VIM_q, Wet$VIM_q)                                        
wilcox.test(Dry$VIM_q, Wet$VIM_q, alternative="greater")
wilcox.test(Dry$VIM_q, Wet$VIM_q, alternative="less")  
###########INT######################
shapiro.test(Dry$INT_q) #0.008578
shapiro.test(Wet$INT_q) #0.00512
#3.1.2 Two Sample t-test (unpaired)
t.test(Dry$INT_q, Wet$INT_q)
#one tail test
t.test(Dry$INT_q, Wet$INT_q, paired = TRUE, alternative="greater")
t.test(Dry$INT_q, Wet$INT_q, paired = TRUE, alternative="less")
###########TETQ######################
shapiro.test(Dry$TETQ_q) #0.4627
shapiro.test(Wet$TETQ_q) #0.0129
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$TETQ_q, Wet$TETQ_q)                                        
wilcox.test(Dry$TETQ_q, Wet$TETQ_q, alternative="greater")
wilcox.test(Dry$TETQ_q, Wet$TETQ_q, alternative="less") 
###########SUL#######################
shapiro.test(Dry$SUL_q) #2.641e-05
shapiro.test(Wet$SUL_q) #0.0005379
#3.1.2 Two Sample t-test (unpaired)
t.test(Dry$SUL_q, Wet$SUL_q)
#one tail test
t.test(Dry$SUL_q, Wet$SUL_q, paired = TRUE, alternative="greater")
t.test(Dry$SUL_q, Wet$SUL_q, paired = TRUE, alternative="less")
###########MCR######################
#Para
#3.1.1 MAXIMUM LIKELIHOOD ESTIMATION (unpaired) (see topic 9.3, helsel)/ dist= c("gaussian","lognormal")
Qseason$MCR_q <- as.numeric(Qseason$MCR_q)
Qseason$MCR_cen <- as.logical(Qseason$MCR_cen)
Qseason$Season <- as.factor(Qseason$Season)
tcemle <- cenmle(Qseason$MCR_q, Qseason$MCR_cen, Qseason$Season, dist="lognormal")
tcemle
shapiro.test(Dry$MCR_q) #0.1392
shapiro.test(Wet$MCR_q) #0.01236
#3.1.2 Two Sample t-test (unpaired)
t.test(Dry$MCR_q, Wet$MCR_q)
#one tail test
t.test(Dry$MCR_q, Wet$MCR_q, paired = TRUE, alternative="greater")
t.test(Dry$MCR_q, Wet$MCR_q, paired = TRUE, alternative="less")

###################new qpcr analysis by phage and arg in qpcr#######################

library(readr)
Seasondiff_qpcr <- read_csv("Seasondiff_qpcr.csv", 
                            col_types = cols(crAssphage = col_number(), 
                                             INT_q = col_number(), SUL_q = col_number(), 
                                             TEM_q = col_number(), NDM_q = col_number(), 
                                             VIM_q = col_number(), TETQ_q = col_number(), 
                                             MCR_q = col_number()))
head(Seasondiff_qpcr)
View(Seasondiff_qpcr)

Dry <- Seasondiff_qpcr[1:9,]
Wet <- Seasondiff_qpcr[10:18,]

Seasondiff_qpcr$Season <- as.factor(Seasondiff_qpcr$Season)
sum.stat <- cenfit(Seasondiff_qpcr$MCR_q, Seasondiff_qpcr$MCR_cen, Seasondiff_qpcr$Season) 
sum.stat


########Phage###########
summary(Dry$crAssphage)
sd(Dry$crAssphage)
shapiro.test(Dry$crAssphage)
summary(Wet$crAssphage)
sd(Wet$crAssphage)
shapiro.test(Wet$crAssphage)

#####Phage########
shapiro.test(Dry$crAssphage) #0.1338
shapiro.test(Wet$crAssphage) #0.01284
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$crAssphage, Wet$crAssphage)                                        #0.6665
wilcox.test(Dry$crAssphage, Wet$crAssphage, alternative="greater")                 #0.3332
wilcox.test(Dry$crAssphage, Wet$crAssphage, alternative="less")                    #0.6976

####################################################################################
#Ht-qPCR

library(readr)
Hseason <- read_csv("Hseason.csv", col_types = cols(INT_h = col_number(), 
                                                    SUL_h = col_number(), TEM_h = col_number(), 
                                                    NDM_h = col_number(), VIM_h = col_number(), 
                                                    TETQ_h = col_number(), MCR_h = col_number()))


head(Hseason)

Dry <- Hseason[1:9,]
Wet <- Hseason[10:18,]

###########NDM#############
shapiro.test(Dry$NDM_h) #0.09216
shapiro.test(Wet$NDM_h) #0.8672
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$NDM_h, Wet$NDM_h)                                        #0.6665
wilcox.test(Dry$NDM_h, Wet$NDM_h, alternative="greater")                 #0.3332
wilcox.test(Dry$NDM_h, Wet$NDM_h, alternative="less")                    #0.6976
###########TEM######################
shapiro.test(Dry$TEM_h) #0.001507
shapiro.test(Wet$TEM_h) #0.3932
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$TEM_h, Wet$TEM_h)                                        
wilcox.test(Dry$TEM_h, Wet$TEM_h, alternative="greater")
wilcox.test(Dry$TEM_h, Wet$TEM_h, alternative="less")  
###########VIM######################
shapiro.test(Dry$VIM_h) #0.3627
shapiro.test(Wet$VIM_h) #0.399
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$VIM_h, Wet$VIM_h)                                        
wilcox.test(Dry$VIM_h, Wet$VIM_h, alternative="greater")
wilcox.test(Dry$VIM_h, Wet$VIM_h, alternative="less")  
###########INT######################
shapiro.test(Dry$INT_h) #0.3292
shapiro.test(Wet$INT_h) #0.5525
#3.1.2 Two Sample t-test (unpaired)
t.test(Dry$INT_h, Wet$INT_h)
#one tail test
t.test(Dry$INT_h, Wet$INT_h, paired = TRUE, alternative="greater")
t.test(Dry$INT_h, Wet$INT_h, paired = TRUE, alternative="less")
###########TETQ######################
shapiro.test(Dry$TETQ_h) #0.1453
shapiro.test(Wet$TETQ_h) #0.9961
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$TETQ_h, Wet$TETQ_h)                                        
wilcox.test(Dry$TETQ_h, Wet$TETQ_h, alternative="greater")
wilcox.test(Dry$TETQ_h, Wet$TETQ_h, alternative="less") 
###########SUL#######################
shapiro.test(Dry$SUL_h) #0.02573
shapiro.test(Wet$SUL_h) #0.1207
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$SUL_h, Wet$SUL_h)                                        
wilcox.test(Dry$SUL_h, Wet$SUL_h, alternative="greater")
wilcox.test(Dry$SUL_h, Wet$SUL_h, alternative="less") 
###########MCR######################
#Para
#3.1.1 MAXIMUM LIKELIHOOD ESTIMATION (unpaired) (see topic 9.3, helsel)/ dist= c("gaussian","lognormal")
Hseason$MCR_h <- as.numeric(Hseason$MCR_h)
Hseason$MCR_cen <- as.logical(Hseason$MCR_cen)
Hseason$Season <- as.factor(Hseason$Season)

shapiro.test(Dry$MCR_h) #0.0001948
shapiro.test(Wet$MCR_h) #0.001054
#3.1.2 Two Sample t-test (unpaired)
t.test(Dry$MCR_h, Wet$MCR_h)
#one tail test
t.test(Dry$MCR_h, Wet$MCR_h, paired = TRUE, alternative="greater")
t.test(Dry$MCR_h, Wet$MCR_h, paired = TRUE, alternative="less")


##########Site sum qPCR#############
library(readr)
library(readr)
Season_sumdiff_q <- read_csv("Season_sumdiff_q.csv", 
                             col_types = cols(crAssphage = col_number(), 
                                              INT_q = col_number(), SUL_q = col_number(), 
                                              TEM_q = col_number(), NDM_q = col_number(), 
                                              VIM_q = col_number(), TETQ_q = col_number(), 
                                              MCR_q = col_number()))
View(Season_sumdiff_q)
head(Season_sumdiff_q)

WW <- qPCR_sitesumdiff[15:18, ]
AE <- qPCR_sitesumdiff[1:2, ]
CL2 <- qPCR_sitesumdiff[3:6, ]
DW <- qPCR_sitesumdiff[7:10, ]
UP <- qPCR_sitesumdiff[11:14, ]

library(readr)
dwup <- read_csv("dwup.csv", col_types = cols(VIM = col_number()))
View(dwup)

DW <- dwup[17:20, ]
UP <- dwup[21:24, ]
head(DW)
summary(DW$VIM)
sd(DW$VIM)
shapiro.test(DW$VIM)
summary(UP$VIM)
sd(UP$VIM)
shapiro.test(UP$VIM)

######site sum diff qpcr ######
dwup$Site <- as.factor(dwup$Site)
# 4.1.1 ANOVA
one.way.anova_event <- aov(dwup$VIM ~ dwup$Site, data = dwup)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

#######WW############
summary(WW$TEM_q)
sd(WW$TEM_q)
shapiro.test(WW$TEM_q)

summary(WW$VIM_q)
sd(WW$VIM_q)
shapiro.test(WW$VIM_q)

summary(WW$INT_q)
sd(WW$INT_q)
shapiro.test(WW$INT_q)

summary(WW$MCR_q) #DetectedOnly
sd(WW$MCR_q)
shapiro.test(WW$MCR_q)

summary(WW$SUL_q)
sd(WW$SUL)
shapiro.test(WW$SUL_q)

summary(WW$NDM_q)
sd(WW$NDM_q)
shapiro.test(WW$NDM_q)

summary(WW$TETQ_q)
sd(WW$TETQ_q)
shapiro.test(WW$TETQ_q)

summary(WW$crAssphage)
sd(WW$crAssphage)
shapiro.test(WW$crAssphage)
#########ARG_AE###########
summary(AE$NDM_q)
sd(AE$NDM_q)


summary(AE$TEM_q)
sd(AE$TEM_q)


summary(AE$VIM_q)
sd(AE$VIM_q)


summary(AE$INT_q)
sd(AE$INT_q)


summary(AE$MCR_q) #Non-detected
sd(AE$MCR_q)


summary(AE$SUL_q)
sd(AE$SUL_q)


summary(AE$TETQ_q)
sd(AE$TETQ_q)


summary(AE$crAssphage)
sd(AE$crAssphage)


#########ARG_CL2###########
summary(CL2$NDM_q)
sd(CL2$NDM_q)
shapiro.test(CL2$NDM_q)

summary(CL2$TEM_q)
sd(CL2$TEM_q)
shapiro.test(CL2$TEM_q)

summary(CL2$VIM_q)
sd(CL2$VIM_q)
shapiro.test(CL2$VIM_q)

summary(CL2$INT_q)
sd(CL2$INT_q)
shapiro.test(CL2$INT_q)

summary(CL2$MCR_q) #Non-detected
sd(CL2$MCR_q)
shapiro.test(CL2$MCR_q)

summary(CL2$SUL_q)
sd(CL2$SUL_q)
shapiro.test(CL2$SUL_q)

summary(CL2$TETQ_q)
sd(CL2$TETQ_q)
shapiro.test(CL2$TETQ_q)

summary(CL2$crAssphage)
sd(CL2$crAssphage)
shapiro.test(CL2$crAssphage)

#########ARG_DW###########
summary(DW$NDM_q)
sd(DW$NDM_q)
shapiro.test(DW$NDM_q)

summary(DW$TEM_q)
sd(DW$TEM_q)
shapiro.test(DW$TEM_q)

summary(DW$VIM_q)
sd(DW$VIM_q)
shapiro.test(DW$VIM_q)

summary(DW$INT_q)
sd(DW$INT_q)
shapiro.test(DW$INT_q)

summary(DW$MCR_q) #Non-detected
sd(DW$MCR_q)
shapiro.test(DW$MCR_q)

summary(DW$SUL_q)
sd(DW$SUL_q)
shapiro.test(DW$SUL_q)

summary(DW$TETQ_q)
sd(DW$TETQ_q)
shapiro.test(DW$TETQ_q)

summary(DW$crAssphage)
sd(DW$crAssphage)
shapiro.test(DW$crAssphage)

#########ARG_UP###########
summary(UP$NDM_q)
sd(UP$NDM_q)
shapiro.test(UP$NDM_q)

summary(UP$TEM_q)
sd(UP$TEM_q)
shapiro.test(UP$TEM_q)

summary(UP$VIM_q)
sd(UP$VIM_q)
shapiro.test(UP$VIM_q)

summary(UP$INT_q)
sd(UP$INT_q)
shapiro.test(UP$INT_q)

summary(UP$MCR_q) #Non-detected
sd(UP$MCR_q)
shapiro.test(UP$MCR_q)

summary(UP$SUL_q)
sd(UP$SUL_q)
shapiro.test(UP$SUL_q)

summary(UP$TETQ_q)
sd(UP$TETQ_q)
shapiro.test(UP$TETQ_q)

summary(UP$crAssphage)
sd(UP$crAssphage)
shapiro.test(UP$crAssphage)

######site sum diff qpcr ######
Season_sumdiff_q$Site <- as.factor(Season_sumdiff_q$Site)
# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$NDM_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$TEM_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$VIM_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$INT_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$MCR_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$SUL_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$TETQ_q ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(Season_sumdiff_q$crAssphage ~ Season_sumdiff_q$Site, data = Season_sumdiff_q)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

######sum site diff ht########
library(readr)
HT_sitesum <- read_csv("HT_sitesum.csv", 
                       col_types = cols(INT_h = col_number(), 
                                        SUL_h = col_number(), TEM_h = col_number(), 
                                        NDM_h = col_number(), VIM_h = col_number(), 
                                        TETQ_h = col_number(), MCR_h = col_number()))
View(HT_sitesum)
head(HT_sitesum)

WW <- HT_sitesum[1:4, ]
AE <- HT_sitesum[5:6, ]
CL2 <- HT_sitesum[7:10, ]
DW <- HT_sitesum[11:14, ]
UP <- HT_sitesum[15:18, ]
#######WW############
summary(WW$Phage)
sd(WW$Phage)
shapiro.test(WW$Phage)

summary(AE$Phage)
sd(AE$Phage)
shapiro.test(AE$Phage)

summary(CL2$Phage)
sd(CL2$Phage)
shapiro.test(CL2$Phage)

summary(DW$Phage)
sd(DW$Phage)
shapiro.test(DW$Phage)

summary(UP$Phage)
sd(UP$Phage)
shapiro.test(UP$Phage)

summary(WW$TEM_h)
sd(WW$TEM_h)
shapiro.test(WW$TEM_h)

summary(WW$VIM_h)
sd(WW$VIM_h)
shapiro.test(WW$VIM_h)

summary(WW$INT_h)
sd(WW$INT_h)
shapiro.test(WW$INT_h)

summary(WW$MCR_h)
sd(WW$MCR_h)
shapiro.test(WW$MCR_h)

summary(WW$SUL_h)
sd(WW$SUL_h)
shapiro.test(WW$SUL_h)

summary(WW$NDM_h)
sd(WW$NDM_h)
shapiro.test(WW$NDM_h)

summary(WW$TETQ_h)
sd(WW$TETQ_h)
shapiro.test(WW$TETQ_h)

#########ARG_AE###########
summary(AE$NDM_h)
sd(AE$NDM_h)


summary(AE$TEM_h)
sd(AE$TEM_h)


summary(AE$VIM_h)
sd(AE$VIM_h)


summary(AE$INT_h)
sd(AE$INT_h)


summary(AE$MCR_h) #Non-detected
sd(AE$MCR_h)


summary(AE$SUL_h)
sd(AE$SUL_h)

summary(AE$TETQ_h)
sd(AE$TETQ_h)


#########ARG_CL2###########
summary(CL2$NDM_h)
sd(CL2$NDM_h)
shapiro.test(CL2$NDM_h)

summary(CL2$TEM_h)
sd(CL2$TEM_h)
shapiro.test(CL2$TEM_h)

summary(CL2$VIM_h)
sd(CL2$VIM_h)
shapiro.test(CL2$VIM_h)

summary(CL2$INT_h)
sd(CL2$INT_h)
shapiro.test(CL2$INT_h)

summary(CL2$MCR_h) #Non-detected
sd(CL2$MCR_h)
shapiro.test(CL2$MCR_h)

summary(CL2$SUL_h)
sd(CL2$SUL_h)
shapiro.test(CL2$SUL_h)

summary(CL2$TETQ_h)
sd(CL2$TETQ_h)
shapiro.test(CL2$TETQ_h)

#########ARG_DW###########
summary(DW$NDM_h)
sd(DW$NDM_h)
shapiro.test(DW$NDM_h)

summary(DW$TEM_h)
sd(DW$TEM_h)
shapiro.test(DW$TEM_h)

summary(DW$VIM_h)
sd(DW$VIM_h)
shapiro.test(DW$VIM_h)

summary(DW$INT_h)
sd(DW$INT_h)
shapiro.test(DW$INT_h)

summary(DW$MCR_h) 
sd(DW$MCR_h)
shapiro.test(DW$MCR_h)

summary(DW$SUL_h)
sd(DW$SUL_h)
shapiro.test(DW$SUL_h)

summary(DW$TETQ_h)
sd(DW$TETQ_h)
shapiro.test(DW$TETQ_h)

#########ARG_UP###########
summary(UP$NDM_h)
sd(UP$NDM_h)
shapiro.test(UP$NDM_h)

summary(UP$TEM_h)
sd(UP$TEM_h)
shapiro.test(UP$TEM_h)

summary(UP$VIM_h)
sd(UP$VIM_h)
shapiro.test(UP$VIM_h)

summary(UP$INT_h)
sd(UP$INT_h)
shapiro.test(UP$INT_h)

summary(UP$MCR_h)
sd(UP$MCR_h)
shapiro.test(UP$MCR_h)

summary(UP$SUL_h)
sd(UP$SUL_h)
shapiro.test(UP$SUL_h)

summary(UP$TETQ_h)
sd(UP$TETQ_h)
shapiro.test(UP$TETQ_h)

######site sum diff qpcr ######
HT_sitesum$Site <- as.factor(HT_sitesum$Site)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$NDM_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$TEM_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$VIM_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$INT_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$MCR_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$SUL_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$TETQ_h ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)

# 4.1.1 ANOVA
one.way.anova_event <- aov(HT_sitesum$Phage ~ HT_sitesum$Site, data = HT_sitesum)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)


# phage absolute Season_sum qPCR

library(readr)
sum_phage_ab <- read_csv("sum_phage_ab.csv", 
                         col_types = cols(crAssphage = col_number()))
View(sum_phage_ab)
head(sum_phage_ab)

Dry <- sum_phage_ab[1:12,]
Wet <- sum_phage_ab[13:24,]

sum_phage_ab$Season <- as.factor(sum_phage_ab$Season)

########Phage###########
summary(Dry$crAssphage)
sd(Dry$crAssphage)
shapiro.test(Dry$crAssphage)
summary(Wet$crAssphage)
sd(Wet$crAssphage)
shapiro.test(Wet$crAssphage)

#####Phage########
shapiro.test(Dry$crAssphage) #0.8036
shapiro.test(Wet$crAssphage) #0.4576
#3.2.2 Wilcoxon rank sum test or Mann-Whitney U Test (unpaired)
wilcox.test(Dry$crAssphage, Wet$crAssphage)                                        #0.2913
wilcox.test(Dry$crAssphage, Wet$crAssphage, alternative="greater")                 #0.1457
wilcox.test(Dry$crAssphage, Wet$crAssphage, alternative="less")                    #0.8672

##########Site sum qPCR#############
library(readr)
sumsite_phageab <- read_csv("sumsite_phageab.csv", 
                            col_types = cols(crAssphage = col_number()))
View(sumsite_phageab)
head(sumsite_phageab)

WW <- sumsite_phageab[1:4, ]
AE <- sumsite_phageab[5:8, ]
CL1 <- sumsite_phageab[9:12, ]
CL2 <- sumsite_phageab[13:16, ]
DW <- sumsite_phageab[17:20, ]
UP <- sumsite_phageab[21:24, ]

summary(WW$crAssphage)
sd(WW$crAssphage)
shapiro.test(WW$crAssphage)

summary(AE$crAssphage)
sd(AE$crAssphage)
shapiro.test(AE$crAssphage)

summary(CL1$crAssphage)
sd(CL1$crAssphage)
shapiro.test(CL1$crAssphage)

summary(CL2$crAssphage)
sd(CL2$crAssphage)
shapiro.test(CL2$crAssphage)

summary(DW$crAssphage)
sd(DW$crAssphage)
shapiro.test(DW$crAssphage)

summary(UP$crAssphage)
sd(UP$crAssphage)
shapiro.test(UP$crAssphage)

# 4.1.1 ANOVA
sumsite_phageab$Site <- as.factor(sumsite_phageab$Site)
one.way.anova_event <- aov(sumsite_phageab$crAssphage ~ sumsite_phageab$Site, data = sumsite_phageab)
# Summary of the analysis
summary(one.way.anova_event)
#Tukey's multiple comparisons test
TukeyHSD(one.way.anova_event)
