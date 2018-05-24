#These are the packages you need for the analysis
``
library(plyr)
library(rentrez)
``
#Read the table from magic-blast in
tab=read.table("Magic-BLAST_output_example",header=F)
#fixing column names
col_names=c("query.acc.","reference.acc.","% identity","not used","not.used",
            "not.used","query.start","query.end","reference.start","reference.end","not.used",
            "not.used","score","query.strand","reference.strand","query.length","BTOP",
            "num.placements","not.used","compartment","left.overhang","right.overhang",
            "mate.reference","mate.ref..start","composite.score")
colnames(tab)=col_names
#remove duplicate hits for eahc sequence
tab$seq=sub("..$","",tab$query.acc.)
tab2=tab[!duplicated(tab[,"seq"]),]
#clean SRA names to aggregate them
tab2$SRA=sub("\\..*","",tab2$seq)
#tab2$length=abs(tab2$query.start-tab2$query.end)
#tab2=tab2[tab2$length>=50,]
#tab2=tab2[tab2$X..identity>=90,]
#count freq of each enzyme for SRA
df <- count(tab2, c('SRA','reference.acc.'))
#save the clean table
# write.csv(df,"tablegene.csv")
#fuction for retrieving SRA MBS size
getSize <- function(ids) {size_mega <- c()
  for(i in ids) {term =  paste0(i, "[ACCN]")
    run = entrez_search(db = "sra", term = term)
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    ## extract range defined by "total bases" and "load_done"
    size = substr(x, start = regexpr("total_bases=", x)[[1]][1] + attr(regexpr("total_bases=", x), "match.length"),stop = regexpr(" load_done", x)[[1]][1])
    ### some clean up
    size = gsub('\"', "", size, fixed = TRUE)
    ## converto to numeric and Mb
    size = as.numeric(size)/1e6
    size_mega <- c(size_mega, size)}
size_mega}
#retrieve the size of SRAs
df$MBS=getSize(df$SRA)
#Noprmalize using the size of each SRA
df$count=(df$freq*mean(df$MBS)/df$MBS)
#Transform to log10
df$trans_count=log10(df$count)
#Divide in tertiles or percentiles
df$group <- as.numeric(cut(df$trans_count, 4))



#Plotting boxplot of the pecentiles
g <- ggplot(df, aes(group, trans_count))
g + geom_boxplot(aes(fill=factor(group)))+ scale_fill_brewer(palette="Set1")+guides(fill=FALSE)+
   labs(x="Female Gut Metagenome",y="log10 of Beta-Glucuronidase Hits")+theme_classic(base_size = 30)+
   geom_point()+scale_x_continuous(breaks=c(1, 2,3,4),labels=c("Low", "Medium","High","Very High"))
