<img src="./Logo.png">

# PYRAMMID (PRedict Metagenome Medication Interaction)

A simple pipeline to predict metagenome medication Interaction.

# Bioinformatics Analysis
# Structural similarity analysis of chemical compounds
## MCS analysis
```
library(ChemmineR)
library(fmcsR)
sdfset <- read.SDFset("compounds.sdf")
tanimoto_coeff <- sapply(cid(sdfset), function(x) fmcsBatch(sdfset[x], sdfset, matching.mode="aromatic",au=1, bu=1,numParallel = 6)[,"Tanimoto_Coefficient"])
write.csv(tanimoto_coeff, "fcms-Tanimoto_Coefficient.csv")
```
## Network analysis visualization
```
library(reshape)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
df=read.csv("fcms-Tanimoto_Coefficient.csv",header=T,row.names = 1)
df=as.matrix(df)
n=rownames(df)
upperTriangle=upper.tri(df, diag=F)
cor.upperTriangle=df
cor.upperTriangle[!upperTriangle]=NA
cor_melted<-na.omit(melt(cor.upperTriangle))
colnames(cor_melted)<-c("X1", "X2", "Tanimoto_Coefficient")
cor_melted=cor_melted[cor_melted$Tanimoto_Coefficient >=.65,]
cor.graph <- as_tbl_graph(cor_melted, directed = FALSE)
cor.group <- data_frame(name = rownames(df))
cor.graph <- cor.graph %>%
   activate(nodes) %>%
   left_join(cor.group, by = "name") %>%
   rename(label = name)
cor.graph %>%
   activate(nodes) %>%
   mutate(Clusters = as.factor(group_infomap())) %>% 
   ggraph(layout = "graphopt") + 
   geom_edge_arc(aes(width = Tanimoto_Coefficient), colour="lightblue") +
   geom_node_point(aes(colour = Clusters), size = 7) +
   geom_node_text(aes(label = label), size=5,repel = TRUE)+
   labs(edge_width ="Tanimoto Coefficient")+
   theme_graph()
```

# Preparation of human gut metagenome data sets
## NCBI SRA Toolkit
```
fastq-dump SRA --read-filter pass --skip-technical --clip --minReadLen 50 --readids --split-3 --outdir SRA
```
## Fastq-join
```
fastq-join SRA_pass_1.fastq SRA_pass_2.fastq -o SRA_%.fastq
```
## PRINSEQ++
```
prinseq++ -fastq SRA -min_qual_mean 20 -ns_max_n 0 -derep -trim_qual_right=20 -lc_entropy -min_len 50 -threads 36 -out_format 1 -out_name SRA
```
## MetaPhyler
```
runMetaphyler.pl SRA.fasta blastn SRA 1
```
## DIAMOND
```
diamond makedb --in reference.fasta -d reference
diamond blastx -d diamond_ref/reference -q SRA -o SRA_matches.m8 -p 36
```

## DIAMOND downstream analysis
```
#Reading the DIAMOND output table
library(plyr)
tab=read.table("diamond_results",header=F)
#Adding column names
colnames(tab)[1]="query.acc."
colnames(tab)[2]="reference.acc."
#Remove duplicate hits for each sequence
tab$reference.acc.=sub("[0-9]","",tab$reference.acc.)
tab$reference.acc.=sub("[0-9]","",tab$reference.acc.)
tab$seq=sub("..$","",tab$query.acc.)
tab2=tab[!duplicated(tab[,"seq"]),]
#Clean SRA names to aggregate them
tab2$SRA=sub("\\..*","",tab2$seq)
df <- count(tab3, c('SRA','reference.acc.'))
#Reading in the MetaPhyler results
taxa=read.table("result_taxa.txt",header=F)
colnames(taxa)=c("SRA","phyla")
almost_final=merge(df, taxa, all.x = TRUE)
#Reading in the metadata for the SRA
metadata=read.csv("metadata.csv",header=T)
final=merge(almost_final, metadata, all.x = TRUE)
#Normalizing the data to the metagenome size and transforming (log2)
final$norm_count=(final$freq*1000)/final$phyla
#Writing everything to a table
write.csv(final2,"results_clean_diamond.csv")
```
