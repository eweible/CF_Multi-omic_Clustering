###############################
##
## Project: Brandie CF data project
##
## Purpose: Cleaning each time point and saving for analysis
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2023-01-26
##
## ---------------------------
## Notes:
##    Do each CF time point with both methods
##
## --------

# Helpful functions
`%nin%` <- Negate(`%in%`)
library(tidyverse); library(magrittr); library(plotly)
library(compositions); library(VIM)
library(MAI); library(ClusterR); library(RVA)
library(ggrepel)

source("C:/Users/ezraw/OneDrive/Documents/Anschutz/RA/Promise_Microbiome_CharlieCarpenter/Promise_Data/Charlie_Files/Carpenter dissertation application/CF_for_Brandie/CF_for_Brandie/SimulationFunctions_Ezra.R")

## Data Read In ----

mb_dat <- read.table("C:/Users/ezraw/OneDrive/Documents/Anschutz/RA/Promise_Microbiome_CharlieCarpenter/Promise_Data/Charlie_Files/Carpenter dissertation application/CF_for_Brandie/CF_for_Brandie/Data/Raw/CF_Infant/seqDataforKayla.txt",
                     header = T, sep = ',')
prot_dat <- read.csv("C:/Users/ezraw/OneDrive/Documents/Anschutz/RA/Promise_Microbiome_CharlieCarpenter/Promise_Data/Charlie_Files/Carpenter dissertation application/CF_for_Brandie/CF_for_Brandie/Data/Raw/CF_Infant/SomaDataforKayla.txt",
                       header = T)

# TimePoint 1 -------------------------------------------------------------

# * MB --------------------------------------------------------------------

mb_dat %>% filter(time == 1) %$% table(sid)

mb_tax1 <- mb_dat %>% add_count(sid_first) %>% 
  filter(time == 1, n > 1, !is.na(seq_count)) %>%
  dplyr::select(sid, taxa, seq_count) %>% 
  pivot_wider(names_from = taxa, values_from = seq_count,
              values_fn = max)

## Prevalence
prev1 <- apply(mb_tax1, 2, function(x) sum(x>0)/length(x))

## Relative abundance
tots1 <- rowSums(mb_tax1[,-1]) # seq depth
ra1 <- sweep(mb_tax1[,-1], 1, tots1, FUN = '/')

## keep microbs that make up at least 0.1% abundance 
## in at least 1 child
## adding TRUE for sid column
ra.keep1 <- c(T, apply(ra1, 2, function(x) any(x > 0.001)) )

## Must be present in 10% of kids (3 kids)
mb_tax1F <- mb_tax1[, prev1 >= 0.10 & ra.keep1]
mb_tax1F %<>% filter(sid != 49)

mb_clin1 <- mb_dat %>% 
  filter(time == 1, sid %in% mb_tax1F$sid) %>% 
  # , total > 5000) %>% 
  dplyr::select(sid, gender, genotypes1) %>% 
  # mutate(sid = paste0("X", sid)) %>% 
  distinct()

tots1F <- rowSums(mb_tax1F[,-1])

mb_clr1F <- sweep(mb_tax1F[,-1], 1, 1/tots1F, "+") %>% #Taking all of the data excluding the sids and sum up the relative abundance
  apply(1, clr) %>% t   #Normalization piece (centralized ratio)

# * Prots -----------------------------------------------------------------

prots1_raw <- prot_dat %>% 
  filter(time == 1, sid %in% mb_tax1F$sid) %>% 
  dplyr::select(sid, Target, meas) %>% 
  pivot_wider(names_from = Target, values_from = meas,
              values_fn = max)

p_miss1 <- apply(prots1_raw, 2, function(x) sum(is.na(x))/length(x))
prots1_raw <- prots1_raw[, p_miss1 < 0.2]
prots1_l2 <- apply(prots1_raw[,-1], 2, log2)

# * Outliers ----------------------------------------------------------------

### ### ### ### Microbes ### ### ### ### 
mbPCA <- mb_clr1F %>% t %>% cov %>% eigen

as.data.frame(mbPCA$vectors) %>% plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3)

### ### ### ### Protiens ### ### ### ###

protPCA <- prots1_l2 %>% t %>% cov %>% eigen

as.data.frame(protPCA$vectors) %>% plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3)

## Outlier: sid 49
# prots1_raw[which.min(protPCA$vectors[,3]),1]

CF_T1 <- list(prots1_l2=prots1_l2, mb1_clr = mb_clr1F)

# TimePoint 2 -------------------------------------------------------------


# * MB --------------------------------------------------------------------

mb_dat %>% filter(time == 2) %$% table(sid)

mb_tax2 <- mb_dat %>% add_count(sid_first) %>% 
  filter(time == 2, n > 1, !is.na(seq_count)) %>%
  dplyr::select(sid, taxa, seq_count) %>% 
  pivot_wider(names_from = taxa, values_from = seq_count,
              values_fn = max)

## Prevalence
prev2 <- apply(mb_tax2, 2, function(x) sum(x>0)/length(x))

## Relative abundance
tots2 <- rowSums(mb_tax2[,-1]) # seq depth
ra2 <- sweep(mb_tax2[,-1], 1, tots2, FUN = '/')

## keep microbs that make up at least 0.1% abundance 
## in at least 1 child
## adding TRUE for sid column
ra.keep2 <- c(T, apply(ra2, 2, function(x) any(x > 0.001)) )

## Must be present in 10% of kids (3 kids)
mb_tax2F <- mb_tax2[, prev2 >= 0.10 & ra.keep2]
mb_tax2F %<>% filter(sid != 49)

mb_clin2 <- mb_dat %>% 
  filter(time == 2, sid %in% mb_tax2F$sid) %>% 
  # , total > 5000) %>% 
  dplyr::select(sid, gender, genotypes1) %>% 
  # mutate(sid = paste0("X", sid)) %>% 
  distinct()

tots2F <- rowSums(mb_tax2F[,-1])

mb_clr2F <- sweep(mb_tax2F[,-1], 1, 1/tots2F, "+") %>% 
  apply(1, clr) %>% t

# * Prots -----------------------------------------------------------------

prots2_raw <- prot_dat %>% 
  filter(time == 2, sid %in% mb_tax2F$sid) %>% 
  dplyr::select(sid, Target, meas) %>% 
  pivot_wider(names_from = Target, values_from = meas,
              values_fn = max)

p_miss2 <- apply(prots2_raw, 2, function(x) sum(is.na(x))/length(x))
prots2_raw <- prots2_raw[, p_miss2 < 0.2]
prots2_l2 <- apply(prots2_raw[,-1], 2, log2)

# * Outliers ----------------------------------------------------------------

### ### ### ### Microbes ### ### ### ### 
mbPCA <- mb_clr2F %>% t %>% cov %>% eigen

as.data.frame(mbPCA$vectors) %>% plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3)

### ### ### ### Protiens ### ### ### ###

protPCA <- prots2_l2 %>% t %>% cov %>% eigen

as.data.frame(protPCA$vectors) %>% plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3)

## Outlier: sid 49 (again)
# prots2_raw[which.min(protPCA$vectors[,3]),1]

CF_T2 <- list(prots2_l2=prots2_l2, mb2_clr = mb_clr2F)

# TimePoint 3 -------------------------------------------------------------


# * MB --------------------------------------------------------------------

mb_dat %>% filter(time == 3) %$% table(sid)

mb_tax3 <- mb_dat %>% add_count(sid_first) %>% 
  filter(time == 3, n > 1, !is.na(seq_count)) %>%
  dplyr::select(sid, taxa, seq_count) %>% 
  pivot_wider(names_from = taxa, values_from = seq_count,
              values_fn = max)

## Prevalence
prev3 <- apply(mb_tax3, 2, function(x) sum(x>0)/length(x))

## Relative abundance
tots3 <- rowSums(mb_tax3[,-1]) # seq depth
ra3 <- sweep(mb_tax3[,-1], 1, tots3, FUN = '/')

## keep microbs that make up at least 0.1% abundance 
## in at least 1 child
## adding TRUE for sid column
ra.keep3 <- c(T, apply(ra3, 2, function(x) any(x > 0.001)) )

## Must be present in 10% of kids (3 kids)
mb_tax3F <- mb_tax3[, prev3 >= 0.10 & ra.keep3]
mb_tax3F %<>% filter(sid != 55)

mb_clin3 <- mb_dat %>% 
  filter(time == 3, sid %in% mb_tax3F$sid) %>% 
  # , total > 5000) %>% 
  dplyr::select(sid, gender, genotypes1) %>% 
  # mutate(sid = paste0("X", sid)) %>% 
  distinct()


# * Prots -----------------------------------------------------------------

prots3_raw <- prot_dat %>% 
  filter(time == 3, sid %in% mb_tax3F$sid) %>% 
  dplyr::select(sid, Target, meas) %>% 
  pivot_wider(names_from = Target, values_from = meas,
              values_fn = max)

p_miss3 <- apply(prots3_raw, 2, function(x) sum(is.na(x))/length(x))
prots3_raw <- prots3_raw[, p_miss3 < 0.2]
prots3_l2 <- apply(prots3_raw[,-1], 2, log2)

### ### ### ### ### ###

## Removing subjects that aren't present in prots
mb_clr3F <- mb_tax3F %>% 
  filter(sid %in% prots3_raw$sid)

tots3F <- rowSums(mb_clr3F[,-1])

mb_clr3F <- sweep(mb_clr3F[,-1], 1, 1/tots3F, "+") %>% 
  apply(1, clr) %>% t

# * Outliers ----------------------------------------------------------------

### ### ### ### Microbes ### ### ### ### 
mbPCA <- mb_clr3F %>% t %>% cov %>% eigen

as.data.frame(mbPCA$vectors) %>% plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3)

### ### ### ### Protiens ### ### ### ###

protPCA <- prots3_l2 %>% t %>% cov %>% eigen

as.data.frame(protPCA$vectors) %>% plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3)

## Outlier: sid 55 
# prots3_raw[which.max(protPCA$vectors[,3]),1]

CF_T3 <- list(prots3_l2=prots3_l2, mb3_clr = mb_clr3F)

# Saving ------------------------------------------------------------------

CFclean <- list(CF_T1=CF_T1, CF_T2=CF_T2, CF_T3=CF_T3)

#save(CFclean, file="C:/Users/ezraw/OneDrive/Documents/Anschutz/RA/Promise_Microbiome_CharlieCarpenter/Promise_Data/Charlie_Files/Carpenter dissertation application/CF_for_Brandie/CF_for_Brandie/RealData/CleanDat/CF/CFcleanAllTimes_Ezra.RData")

write.csv(CFclean, file = "CF_clean_allTimes.csv")


