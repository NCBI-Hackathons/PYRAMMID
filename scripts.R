Bioinformatics Analysis



### Unix command line

##NCBI SRA Toolkit
fastq-dump SRA --read-filter pass --skip-technical --clip --minReadLen 50 --readids --split-3 --outdir SRA

##Fastq-join
fastq-join SRA_pass_1.fastq SRA_pass_2.fastq -o SRA_%.fastq

##PRINSEQ++
prinseq++ -fastq SRA -min_qual_mean 20 -ns_max_n 0 -derep -trim_qual_right=20 -lc_entropy -min_len 50 -threads 36 -out_format 1 -out_name SRA

##MetaPhyler
runMetaphyler.pl SRA.fasta blastn SRA 1

##Magic-BLAST
makeblastdb -in reference.fasta -parse_seqids -dbtype nucl -out reference
magicblast -query SRA  -db reference -outfmt tabular -no_unaligned -reftype transcriptome -num_threads 36 -score 50 >> output.table

##DIAMOND
diamond makedb --in reference.fasta -d reference
diamond blastx -d diamond_ref/reference -q SRA -o SRA_matches.m8 -p 36

###R programming

##Structural similarity analysis
library("ChemmineR")
library(fmcsR)
sdfset <- read.SDFset("compounds.sdf")
sapply(cid(sdfset), function(x) fmcsBatch(sdfset[x], sdfset, au=0, bu=0,numParallel = 4)[,"Tanimoto_Coefficient"]) 

#Network analysis visualization
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
cor_melted=cor_melted[cor_melted$Tanimoto_Coefficient >=.59999,]
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

##Magic-BLAST downstream analysis
#Reading the magic-blast output table
library(plyr)
magic_tab=read.table("output.table",header=F)
#Adding column names
col_names=c("query.acc.","reference.acc.","% identity","not used","not.used",
        "not.used","query.start","query.end","reference.start","reference.end","not.used",
            "not.used","score","query.strand","reference.strand","query.length","BTOP",
            "num.placements","not.used","compartment","left.overhang","right.overhang",
            "mate.reference","mate.ref..start","composite.score")
colnames(tab)=col_names
#Remove duplicate hits for each sequence
magic_tab$reference.acc.=sub("[0-9]","",magic_tab$reference.acc.)
magic_tab$reference.acc.=sub("[0-9]","",magic_tab$reference.acc.)
magic_tab$seq=sub("..$","",magic_tab$query.acc.)
mod_table=magic_tab[!duplicated(magic_tab[,"seq"]),]
#Clean SRA names to aggregate them
mod_table$SRA=sub("\\..*","",mod_table$seq)
mod_table=mod_table[mod_table$`% identity`>=90,]
count_table<- count(mod_table, c('SRA','reference.acc.'))
#Reading in the MetaPhyler results
taxa=read.table("result_taxa.txt",header=F)
colnames(taxa)=c("SRA","phyla")
almost_final=merge(count_table, taxa, all.x = TRUE)
#Reading in the metadata for the SRA
metadata=read.csv("metadata.csv",header=T)
final=merge(almost_final, metadata, all.x = TRUE)
#Normalizing the data to the metagenome size and transforming (log2)
final$norm_count=((final$freq*mean(final$phyla)*10)/final$phyla)
final$trans_norm_count=log2(final$norm_count)
#Writing everything to a table
write.csv(final,"results_clean.csv")

##DIAMOND downstream analysis
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
tab3=tab2[tab2$V3>=90,]
tab3=tab3[tab3$V4>=50,]
df <- count(tab3, c('SRA','reference.acc.'))
#Reading in the MetaPhyler results
taxa=read.table("result_taxa.txt",header=F)
colnames(taxa)=c("SRA","phyla")
almost_final=merge(df, taxa, all.x = TRUE)
#Reading in the metadata for the SRA
metadata=read.csv("metadata.csv",header=T)
final=merge(almost_final, metadata, all.x = TRUE)
#Normalizing the data to the metagenome size and transforming (log2)
final$norm_count=((final$freq*mean(final$phyla)*10)/final$phyla)
final$trans_norm_count=log2(final2$norm_count)
#Writing everything to a table
write.csv(final2,"results_clean_diamond90_50.csv")
