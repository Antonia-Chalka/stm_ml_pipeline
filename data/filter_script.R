library(tidyverse)
set.seed(100)


args = commandArgs(trailingOnly=TRUE)
# test if there is at least 4 argumentS: if not, return an error
if (length(args)!=7) {
  length(args)
  stop("Exactly 7 arguments must be supplied (input file).n", call.=FALSE)
} 

########################### Input Data  ###########################
snpdist_base <- read.table(args[1],sep="\t")
metadata <- read.csv(args[2]) %>% 
  rename(Filename=args[3],Source.Host=args[4],Collection.Year=args[5],Region=args[6]) %>% 
  select(Region,Source.Host,Collection.Year,Filename)
########################### Find clonal copies ###########################

snpdist_base_ident <- snpdist_base %>% 
  filter(!duplicated(paste0(pmax(V1, V2), pmin(V1, V2)))) %>% 
  filter(V1!=V2 & ( V1!="Reference" & V2!="Reference")) %>% 
  filter(V3<=strtoi(args[7])) 

if  (nrow(snpdist_base_ident)==0) {
  print("No clusters detected")
  quit(status=0)
}

snpdist_base_ident_metadata <-snpdist_base_ident %>% 
  left_join(metadata, by=c("V1" ="Filename")) %>%
  rename(V1.Region = Region, V1.Collection.Year = Collection.Year, V1.Source.Host = Source.Host) %>%
  left_join(metadata, by=c("V2" ="Filename")) %>%
  rename(V2.Region = Region, V2.Collection.Year = Collection.Year, V2.Source.Host = Source.Host) %>%
  filter(V1.Region == V2.Region & V1.Collection.Year == V2.Collection.Year & V1.Source.Host == V2.Source.Host)


# Get unique clusters
snpdist_base_clust <- snpdist_base_ident_metadata %>% 
  group_by(V1.Region, V1.Collection.Year, V1.Source.Host) %>% 
  tally()

# The big boy loop
for (i in c(1:nrow(snpdist_base_clust))){
  
  # Get assembly list
  pull1 <- snpdist_base_ident_metadata %>% 
    filter(V1.Region == unlist(snpdist_base_clust[i,1]), V1.Collection.Year == unlist(snpdist_base_clust[i,2]), 
           V1.Source.Host == unlist(snpdist_base_clust[i,3])) %>%
    pull(V1)
  pull2 <- snpdist_base_ident_metadata %>% 
    filter(V1.Region == unlist(snpdist_base_clust[i,1]), V1.Collection.Year == unlist(snpdist_base_clust[i,2]), 
           V1.Source.Host == unlist(snpdist_base_clust[i,3])) %>%
    pull(V2)
  
  # Merge and get unique
  merged <- c(pull1,pull2) 
  merged_unique <- unique(unlist(merged))
  
  snpdist_base_clust[i,4] <- length(merged_unique)
  
  # Print out to file
  filename <- paste(unlist(snpdist_base_clust[i,1]) ,unlist(snpdist_base_clust[i,2]), unlist(snpdist_base_clust[i,3]),"list", sep=".")
  write.table(merged_unique, file=filename, quote=FALSE, sep="",row.names = FALSE,col.names = FALSE)
}
########################### Cluster Plots ###########################
base_cluster_plot <- ggplot(snpdist_base_clust,aes(V1.Collection.Year,n, fill=V1.Region)) +
  geom_bar(stat="identity",color="black") +
  labs(title="Clonal Outbreak Clusters", y = "# of Sequences", x = "Year", fill = "State", subtitle="Clonal Outbreaks defined as sequences with <= 10 core SNP difference occuring at the same state, region, and host") +
  facet_wrap(~V1.Source.Host) +
  theme_bw()
ggsave("base_cluster_static.png", plot= base_cluster_plot,device=png(),height=10,width=20)
