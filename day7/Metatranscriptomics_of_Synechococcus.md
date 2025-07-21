---
title: Metatranscriptomics of Synechococcus

---

# Metatranscriptomics of *Synechococcus* Clade II in the Hawaiʻi Diel Sampling (HaDS)


:computer: STAMPS 2025


## :question: Purpose

Identify considerations for the analysis and interpretation of metatranscriptomics data. 

Gain familiarity with tools for examining metatranscriptomics data. 

Analyze patterns of gene expression using raw and normalized data. 

## :sunny: Background


Marine *Synechococcus* (phylum Cyanobacteriota) are a group of common photosynethic bacteria that are most dominant in coastal systems, but are widely distributed across the sunlit ocean. 

As a major primary producer, Synechococcus gene expression is relavant for ocean biogeochemistry, ecological interactions, and microbial evolution. 

ata come from the multi-omic Hawaiʻi Diel Sampling (Tucker et al., in review). We sampled surface seawater from nearshore Kāneʻohe Bay, Oʻahu and the adjacent offshore every 90 min for 48 hours. (more details and data accesssion information [here](https://merenlab.org/data/hads/)).




## :ocean: Metatranscriptomic and metagenomic read recruitment

We are going to examine read recruitment of metatranscriptomic and metagenomic samples from the Hawaiʻi Diel Sampling to a reference genome of marine *Synechococcus* called M16. 


Before we dive in, a brief recap of the analysis steps:

1. Quality control metagenomic/ metatranscriptomic reads 
	1. FastQC- assess quality (Andrews, 2010)
	2. Cutadapt- remove adapter sequences reamining in metatranscriptomes (Martin 2011)
	3. SortMeRNA- computationally remove the ribosomal RNA sequences that match to SILVA databases (Kopylova et al., 2012)
	4. Final quality control of paired end reads with `iu-filter-quality-minoche` (Eren et al., 2013)
2. Use reference genome to make a "contig.db" in anvi'o v 8.0 (Eren et al., 2021)
    1. Genes were called with Prodigal (Hyatt et al., 2010)
    2. Genes were annotated with NCBI’s Clusters of Orthologous Groups (COGs) (Tatusov et al., 2003) and a customized HMM database of KEGG orthologs (KOfams) (Aramaki et al., 2020)
4. Metagenomic and metatranscriptomic reads were recruited to the reference genome with Bowtie 2 (Langmead & Salzberg, 2012)
5. `anvi-profile` function stored coverage and detection statistics and summarized and visualized with `anvi-summarize` and `anvi-interactive`

If you want more details about how these analysis steps were conducted check out the `Preliminary_Analysis.R` file. 


# Examine the metatranscriptomic read recruitment
```
cp -r /opt/shared-2/MTX_MGX_anvio/ ~/

cd MTX_MGX_anvio/

conda activate anvio

tar -zxvf M16_MTX_MGX_MERGED.tar.gz M16_MTX_MGX_MERGED

#lets start an anvio interactive session to visualize
anvi-interactive -c M16-contigs.db -p M16_MTX_MGX_MERGED/PROFILE.db -C DEFAULT -b EVERYTHING --gene-mode

```



:bulb:Once you run the anvi-interactive command you will need to click your anvio link ([here](https://github.com/mblstamps/stamps2025/wiki/Accessing-our-cloud-computers)) to start the session.



1. Click Draw.
2. Green stickies up.
3. In Main, click additional, dendrogram angle, and change 270 to 180 to see sample names. Draw.



We are a looking at a genome wide view of metagenomic (MGX) and metatranscriptomic (MTX) reads that were recruited to a *Synechococcus* genome known as M16 (*Synechococcus* Clade II). Each line in the ring is a gene in the genome, and the height of the bar shows the amount of coverage of that gene in the metatranscriptomic/metagenomic data.

Each ring is a sample. The samples are labelled based on Project (HADS), Date (20210818- August 18, 2018), hour (H1200, noon), sample type (MTX for metatranscriptomes or MGX for metagenomes), and the site of sampling (xHP1). 





5. In Main, look at the layers section. See how the Max coverage is highly variable across samples. Currently we are looking at read recruitment of metagenomic (MGX) and metatranscriptomic (MTX) samples. Let's remove the MGX samples from our view. First color the MTX samples in orange, by typing MTX in the edit attributes for multiple layers and click enter. Some boxes should be checked. Change the color of the second box. Click draw. 
6. Uncheck the boxes. Next remove the MGX sample by typing MGX in the edit attributes for multiple layers and change hieght to 0. 
7. Make sure all boxes are checked. Adjust the max coverage (max) to have a high threshold (100) to have a better view of the data. 
8. Finally, in Main, Display, change the Items Order to Synteny (this is the order of the genes in the genome).


 
## Interpretation
Spend some time examining the data and thinking about the following questions:

How does coverage change across the metatranscriptomic samples? 

How does coverage change across genes within a metatranscriptomic sample? 



## Now let's go to the gene-level

Place your arrow on the the ring. Hold down control and then right click one of the genes in the plot. You will get a drop down menu, then go to "Inspect gene and context". 

Click on the genes to see function. 


What do you notice about the coverage profiles at the gene level? Pay attaention to the y-axis scale across the samples.



## Continue to R to examine the output
Look more closely at the gene/genome level information, here is data output from the metatranscriptomic data we were visualizing (in addition to the metagenomic read recruitment data as well).

The anvi-summarize output also has a lot of data for you 
how to look at the coverage/detection of individual genes and the coverage/detection of genomes use anvi-summarize.

How does the coverage of the genome differ between day and night in the metagenomes and in the metatranscriptomes?

```{r}
library(tidyverse)
library(ggplot2)
genome_coverage_M16=read.delim("MTX_MGX_anvio/mean_coverage.txt")

##reorganize your genome coverage dataframe, group samples by day and night categories, and add in gene functions

genome_coverage_M16=genome_coverage_M16%>%
  pivot_longer(cols = `HADS_20210818_H1200_MGX_xHP1_006`:`HADS_20210819_H0300_MTX_xHP1_026`, 
               names_to = "sample",
               values_to = "coverage")%>%
  separate(sample, sep="_", c("Project", "Date", "Time", "Type", "Site","ID"), remove=FALSE)%>%
  mutate(DayNight = case_when(
    Time %in% c("H1200", "H1330", "H1500") ~ "Day",
    Time %in% c("H0130", "H0300", "H0430") ~ "Night",
    TRUE ~ NA_character_
  ))


ggplot(genome_coverage_M16, aes(x=DayNight, y=coverage, color=Type))+
  geom_point()+
  facet_wrap(~Type)+
  ylab("Genome Coverage")
```

How do individual genes differ between day and night?

```{r}
gene_functions_M16=read.delim("MTX_MGX_anvio/EVERYTHING-gene_calls.txt")
gene_coverage_M16=read.delim("MTX_MGX_anvio/EVERYTHING-gene_coverages.txt")

##reorganize your gene coverage dataframe, group samples by day and night categories, and add in gene functions

gene_coverage_M16=gene_coverage_M16%>%
  pivot_longer(cols = `HADS_20210818_H1200_MGX_xHP1_006`:`HADS_20210819_H0300_MTX_xHP1_026`, 
               names_to = "sample",
               values_to = "coverage")%>%
  separate(sample, sep="_", c("Project", "Date", "Time", "Type", "Site","ID"), remove=FALSE)%>%
  mutate(DayNight = case_when(
    Time %in% c("H1200", "H1330", "H1500") ~ "Day",
    Time %in% c("H0130", "H0300", "H0430") ~ "Night",
    TRUE ~ NA_character_
  ))%>%
  left_join(., gene_functions_M16, by="gene_callers_id")


ggplot(gene_coverage_M16, aes(x=DayNight, y=coverage, color=Type))+
  geom_point()+
  facet_wrap(~Type)+
  ylab("Genome Coverage")
```

Which genes have the highest transcription coverage? 
```{r}
gene_coverage_M16%>%
  subset(Type=="MTX")%>%
  group_by(DayNight, gene_callers_id, KOfam)%>%
  summarise(mean=mean(coverage))%>%
  arrange(desc(mean))%>%
  head(10)
```  


Many metatranscriptomic analysis approaches use read count data instead of coverage. Using anvi-blitz output we can look at read count, detection, and coverage at the gene level. 

Is looking at read count and coverage the same thing?
What does it mean to have a very high read count but low coverage? 

Let's examine the relationship between read count and coverage in the example below for gene caller ID 304 (highlighted in red).


```{r}
blitz_M16=read.delim("MTX_MGX_anvio/M16_MTX_MGX_blitz_OUTPUT.txt")

# add in gene annotations
blitz_M16 <- blitz_M16 %>%
  left_join(., gene_functions_M16, by="gene_callers_id")%>%
  mutate(type = ifelse(grepl("MTX", sample), "MTX", "MGX"))


blitz_M16 %>%
  subset(type == "MTX") %>%
  mutate(
    sample = str_remove(sample, "HADS_"),
    highlight = ifelse(gene_callers_id == 304, "highlight", "normal")
  ) %>%
  ggplot(aes(x = num_mapped_reads, y = mean_cov, color = highlight)) +
  geom_point() +
  facet_wrap(~sample, scales = "free") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  guides(color = "none") 

#closer look at the data from gene highlighted above (gene caller id 304). What do you notice about data provided. 


blitz_M16%>%
  subset(type=="MTX")%>%
  subset(gene_callers_id == "304")%>%
  select("gene_callers_id", "sample","length" , "num_mapped_reads", "mean_cov", "KOfam")%>%
  head()

```

Generally, genes are around ~1,000 bp in bacterial genomes, however there is variation in gene length. Gene 304 is BIG! 

Take a peak at the distribution of gene length for M16 isolate of Synechococcus. 

```{r}
ggplot(blitz_M16, aes(x=length))+geom_histogram()+xlab("gene length (bp)")
```


Length of gene is accounted for when examining about coverage, but not when examining read count.





## :ocean: Examining transcripts of Synechococcus M16 in coastal Kāneʻohe Bay that have been normalized by paired metagenomes.



There are a number of considerations in the analysis of metatranscriptomes, including the normalization of metrantranscriptomes for sequencing depth, gene length, abundance of the organism, and changes in community composition. 

The normalization used on these data provide an estimate of how expression of a gene changes across samples while accounting for changes in abundance and community composition.

To do this, we used the center log ratio (clr) of the gene transcripts from Synechococcus and subtracted the average clr of Synechococccus single-copy gene coverage across samples. 

We obtained single-copy gene coverages from a metagenomic co-assembly and gene transcript coverages from recruiting reads to a reference genome. For more details about how these analysis steps were conducted, check out the `Preliminary_Analysis.R` file.
 

# Let's examine how gene expression changes in *Synechococcus* over a diel cycle

To start, we might want to examine some genes we know should vary in expression over time. 

Circadian clock genes are a good candidate because Cyanobacteria use set of set circadian clock genes (e.g. kaiA, kaiB, and kaiC) to determine their 24-hour rhythmic gene expression. Curious to know more?: https://doi.org/10.1073/pnas.0130099100


```{r}
df=read.csv("MTX_MGX_anvio/CLR_MTX_MGX_M16.csv")

g1=df %>%
  subset(KOfam=="circadian clock protein KaiC")%>%
  mutate(HADS_Universal_Sample_ID = str_remove(HADS_Universal_Sample_ID, "HADS_202108")) %>%
  mutate(HADS_Universal_Sample_ID = str_remove(HADS_Universal_Sample_ID, "_xHP1")) %>%
  ggplot(aes(x=HADS_Universal_Sample_ID, y=CLR_MGX_MTX, color=Site, group=Site))+geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +facet_wrap(~Site, scale="free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                       fill = NA,size = 1))+
  theme(legend.position = "none")+scale_color_manual(values=c(HP1="orange"))+
  geom_smooth(se=FALSE,method = "loess", span = 0.4)+
  ylab("kaiC clr MTX MGX")
g1


```
One of the major environmental drivers that result in daily differences in gene expression is the change in sunlight availability over the 24 hour period. As a photosynthetic organism that utilizes sunlight as its main form of energy, the photosynthetic genes of Synechococcus are likely to vary between night and day. 

Let's take a look.

```{r}
g1=df %>%
  filter(str_detect(KOfam, "photosystem")) %>%
  mutate(HADS_Universal_Sample_ID = str_remove(HADS_Universal_Sample_ID, "HADS_202108")) %>%
  mutate(HADS_Universal_Sample_ID = str_remove(HADS_Universal_Sample_ID, "_xHP1")) %>%
  ggplot(aes(x=HADS_Universal_Sample_ID, y=CLR_MGX_MTX, color=Site, group=Site))+geom_point()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +facet_wrap(~Site+KOfam, scale="free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(color = "black",
                                                                       fill = NA,size = 1))+
  theme(legend.position = "none")+scale_color_manual(values=c(HP1="orange"))+
  geom_smooth(se=FALSE,method = "loess", span = 0.4)+
  ylab("photosystem related genes clr MTX MGX")
g1
ggsave("MTX_MGX_anvio/photosyn_genes_M16.png", height=10, width=20)

```
Odd to think of photosystem genes more active at night relative to the day? 

That's what we thought too, but it turns out the nightly expression of photosynthesis genes is associated with preparation for efficient harvest of engery in the early morning hours. This tracks with the most elevated transcription occuring right before dawn at 4:30am (Ottesen et al., 2015 (doi: 10.1073/pnas.1502883112); Becker et al., 2020 (doi: 10.1038/s41396-020-00793-x)).




Let's see how genes differ in transcription on average between night and day relative to the genome 
abundance.

```{r}
df_day_night <- df %>%
  group_by(KOfam, DayNight, COG20_CATEGORY) %>%
  summarise(mean_CLR_MGX_MTX = mean(CLR_MGX_MTX), n = n()) %>%
  ungroup %>%
  pivot_wider(values_from = c("mean_CLR_MGX_MTX", "n"), names_from = "DayNight") %>%
  mutate(diff = mean_CLR_MGX_MTX_Day - mean_CLR_MGX_MTX_Night)

#### more in night
more_active_during_night = df_day_night %>%
  arrange(diff)%>%
  filter(diff<0)
more_active_during_night

head(more_active_during_night)

# many genes are more highly transcribed at night relative to the metagenome. 



#### more in day
more_active_during_day = df_day_night %>%
  arrange(desc(diff))%>%
  filter(diff>0)
more_active_during_day

head(more_active_during_day)

```

Relative to the genome abundance, no genes are more transcriptionally active during the day than the night. We hypothesize that this pattern is driven in part by a decrease in the abundance of *Synechococcus* between and day, that is supported by both our metagenomic clr estimates and absolute cell abundnaces of *Synechococcus*. 



Recent work also suggests that temperature stress can increase transcriptional activity at night in *Synechococcus* clade II (https://doi.org/10.1016/j.algal.2024.103840). Our samples were taken during a period of abnormally high temperatures. If you have extra time, conduct a search for genes related to heat stress! 






Conduct a t-test to examine the differences in transcription between day and night. Do all genes signficiantly differ between day and night?

```{r}

library(tidyverse)
library(dplyr)
library(broom)
install.packages("broom")
library(broom)
names(df)

results <- df %>%
  group_by(gene_callers_id, KOfam) %>%
  group_split() %>%
  lapply(function(group_data) {
    test <- t.test(CLR_MGX_MTX ~ DayNight, data = group_data)
    tidy_result <- broom::tidy(test)
    tidy_result$gene_callers_id <- unique(group_data$gene_callers_id)
    tidy_result$KOfam <- unique(group_data$KOfam)
    tidy_result
  }) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::select(gene_callers_id, KOfam, estimate, statistic, p.value, p_adj, conf.low, conf.high)

results%>%head(20)

view(results)

```




To examine how *Synechococcus* gene expression changes between day and night more broadly, we can examine average differences in the expression of genes across different functional categories. The COG20 Categories come in handy here.

```{r}
df_day_night %>%
  separate_rows(COG20_CATEGORY, sep = "\\|") %>%
  group_by(COG20_CATEGORY) %>%
  summarise(mean_diff = mean(diff, na.rm=TRUE), n=n(), sd_diff=sd(diff, na.rm=TRUE)) %>%
  arrange(desc(abs(mean_diff))) %>%
  head(20)%>%
  print(n = 20)
```


Finally, which genes are the most highly transcribed on average? 

```{r}
df %>%
  group_by(gene_callers_id, KOfam) %>%
  summarise(mean_CLR_MGX_MTX = mean(CLR_MGX_MTX), n = n()) %>%
  ungroup %>%
  arrange(desc(mean_CLR_MGX_MTX))%>%
  head(40)%>%
  print(n = 40)
```
Similar to the raw mean coverage from the metatranscriptomes we see photosynthesis genes are still some of the most abundant after CLR normalization. Cool cool cool. 



## :pencil: Resources


Recap on common normalization strategies: 

https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

Working with MTX and MGX in anvi'o:

https://merenlab.org/2015/06/10/combining-omics-data/

A review on metatranscriptomics in microbial eukaryotes: 

https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2022.867007/full#B34

Consideration for metatranscriptomics: 

https://www.nature.com/articles/ismej201294





## :question: Questions?

If after the course you have additional questions or are interested in chatting, please reach out :) (stucker@mbl.edu). 

