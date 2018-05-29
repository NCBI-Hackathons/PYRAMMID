#These are the packages you need for the analysis in R
library(plyr)
library(rentrez)

#Read the table from Magic-BLAST (The file is provided in this repository)
tab=read.table("Magic-BLAST_output_example",header=F)

#Here we are changing the column names
col_names=c("query.acc.","reference.acc.","% identity","not used","not.used",
            "not.used","query.start","query.end","reference.start","reference.end","not.used",
            "not.used","score","query.strand","reference.strand","query.length","BTOP",
            "num.placements","not.used","compartment","left.overhang","right.overhang",
            "mate.reference","mate.ref..start","composite.score")
colnames(tab)=col_names

#Since the same sequence will have multiple hits, we will remove duplicate hits for each sequence
tab$seq=sub("..$","",tab$query.acc.)
tab2=tab[!duplicated(tab[,"seq"]),]

#This cleans the SRA names to aggregate them before counting
tab2$SRA=sub("\\..*","",tab2$seq)

#Since some hits are very short (20 bases) and some long (100), it is up to you to decide if you want to remove the short ones
tab2$length=abs(tab2$query.start-tab2$query.end)
tab2=tab2[tab2$length>=50,]

#Also, this to filter the identity percentage for the hits, this is up to you
tab2=tab2[tab2$X..identity>=90,]

#To count how many hits per SRA
df <- count(tab2, c('SRA','reference.acc.'))

#If you want to save the cleaned table
write.csv(df,"tablegene.csv")

#To normalize the counts based on the size of SRA, this is a function for retrieving the size of SRA in MBS (developed by Jose)
getSize <- function(ids) {size_mega <- c()
  for(i in ids) {term =  paste0(i, "[ACCN]")
    run = entrez_search(db = "sra", term = term)
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    size = substr(x, start = regexpr("total_bases=", x)[[1]][1] + 
    attr(regexpr("total_bases=", x), "match.length"),
    stop = regexpr(" load_done", x)[[1]][1])
    size = gsub('\"', "", size, fixed = TRUE)
    size = as.numeric(size)/1e6
    size_mega <- c(size_mega, size)}
size_mega}

#To retrieve the MBS for each SRA in a new column
df$MBS=getSize(df$SRA)

#Normalize using the size of each SRA, this is only one way to do it and normalizes only between samples
df$count=(df$freq*mean(df$MBS)/df$MBS)

#Here I transform the numbers using log10 because the variance is high
df$trans_count=log10(df$count)

#Plotting boxplot of all SRAs
g = ggplot(df, aes(reference.acc., trans_count))
g + geom_boxplot() + guides(fill=FALSE) + labs(x="Female Gut Metagenome",y="log10 of Beta-Glucuronidase Hits") +
theme_classic(base_size = 20) + geom_point()

![](https://github.com/NCBI-Hackathons/PYRAMMID/blob/master/output_example.png)

#Here is a very good tutorial to follow depnding on the data you have and if you have multiple groups http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html
