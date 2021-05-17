# The script mostly revolves around text mining of PubMed data with PubTator and more
# Author: Jagajjit Sahu
# Date created: 3rd Jan 2021
# Last updated: 17th May 2021

library(tidyr)
library(igraph)

# Set the path of parent directory
setwd("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin")

# Merging all the PMID files into one
pmidFile1 <- read.table("pubmed/pmid-Proteome_1_8717.txt", header = F, sep = "\t")
pmidFile2 <- read.table("pubmed/pmid-Proteome_2_9552.txt", header = F, sep = "\t")
pmidFile3 <- read.table("pubmed/pmid-Proteome_3_8945.txt", header = F, sep = "\t")
pmidFile4 <- read.table("pubmed/pmid-Proteome_4_5814.txt", header = F, sep = "\t")
pmidFile <- rbind(pmidFile1,pmidFile2,pmidFile3,pmidFile4)
colnames(pmidFile) <- "PMID"
write.table(pmidFile, file = "pubmed/pmid-Proteome.txt", sep = "\t", row.names = F)


# Merging all the pubmed files into one
pubmedFile1 <- readLines("pubmed/bibtext-Proteome_1_8717.nbib")
pubmedFile2 <- readLines("pubmed/bibtext-Proteome_2_9552.nbib")
pubmedFile3 <- readLines("pubmed/bibtext-Proteome_3_8945.nbib")
pubmedFile4 <- readLines("pubmed/bibtext-Proteome_4_5814.nbib")
pubmedFile <- rbind(pubmedFile1,pubmedFile2,pubmedFile3,pubmedFile4)
fileConn<-file("pubmed/pubmed-Proteome.txt")
writeLines(c(pubmedFile1,pubmedFile2,pubmedFile3,pubmedFile4), fileConn)
#writeLines(pubmedFile, fileConn)
close(fileConn)

library(bibliometrix)
# Year-wise growth from nbib data
##Import pubmed citation data in text format and convert it to dataframe 
M <- convert2df("pubmed/pubmed-Proteome.txt", dbsource = "pubmed", format = "pubmed")
#head(M["TC"])
write.csv(M, file = "pubmed/pubmed-Proteome.csv", row.names = F)
write.table(M, file = "pubmed/pubmed-Proteome.txt", sep = "\t", row.names = F)
##For year-wise graph a file with years and corresponding no (also fractionated) of articles
yearWisePub <- as.data.frame(table(M$PY))
colnames(yearWisePub)[1] <- "Year"
###fractionated is the actual no / 1000
yearWisePub$FreqFractnd <- yearWisePub$Freq/1000
##Plot
Publication_Years <- yearWisePub$Year
Freq <- yearWisePub$Freq
Freq_Fractn <- yearWisePub$FreqFractnd

tiff("pubmed/yearWisePub.tiff",height=1800,width=800, res = "300")
barplot(height=Freq_Fractn, names=Publication_Years, 
        col="#0099f9",
        horiz=T, las=1
)
dev.off()

library(ggplot2)

tiff("pubmed/yearWisePub2.tiff",height=1800,width=1200, res = "300")
ggplot(yearWisePub,aes(x=Freq_Fractn, y=Publication_Years)) +
  geom_col(fill="#0099f9") +
  geom_text(aes(label=Freq), vjust=0, size=2)
coord_flip()
dev.off()


# Function to extract pubtator files in a directory and output a merged Entity info file
xtrctPubTat <- function(inDir = getwd(), outDir = getwd()){
  listFiles <- list.files(inDir)
  for (i in 1:length(listFiles)) {
    inptFile <- listFiles[i]
    pubtatDat <- readLines(paste0(inDir, "/", inptFile))
    fileName <- gsub(".*/|.pubtator", "", inptFile)
    
    titl <- grep("[|]t[|]", pubtatDat, value = T)
    abst <- grep("[|]a[|]", pubtatDat, value = T)
    pmid <- gsub("(.*?)\\|(.*)", "\\1", titl)
    titl <- gsub("(.*)\\|(.*)", "\\2", titl)
    abst <- gsub("(.*)\\|(.*)", "\\2", abst)
    IdTitlAbst.df <- data.frame(PMID=pmid, Title=titl, Abstract=abst)
    write.csv(IdTitlAbst.df, file = paste0(outDir,"/", fileName, "_IdTitlAbst.csv"), row.names = F)
    pubtatDatWithtTitlAbst <- pubtatDat[!grepl("[|]t|[|]a", pubtatDat)]
    pubtatDatWithtTitlAbst <- pubtatDatWithtTitlAbst[pubtatDatWithtTitlAbst != ""]
    df1<- as.data.frame(pubtatDatWithtTitlAbst)
    pubtatDatWithtTitlAbst.df <- separate(df1, pubtatDatWithtTitlAbst, c("PMID","Start","End","Entity","Class","Entity_ID"), sep="\t")
    write.csv(pubtatDatWithtTitlAbst.df, file = paste0(outDir,"/", fileName, "_pubtatEntity.csv"), row.names = F)
    print(fileName)
  }
}

#Run
xtrctPubTat("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubtator", "E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubtator")


# Extraction of the entities and the respective information
pubtatDat <- readLines("pubtator_central_export_8Oct20.pubtator")
titl <- grep("[|]t", pubtatDat, value = T)
abst <- grep("[|]a", pubtatDat, value = T)
pmid <- gsub("(.*?)\\|(.*)", "\\1", titl)
titl <- gsub("(.*)\\|(.*)", "\\2", titl)
abst <- gsub("(.*)\\|(.*)", "\\2", abst)
IdTitlAbst.df <- data.frame(PMID=pmid, Title=titl, Abstract=abst)
write.csv(IdTitlAbst.df, file = "IdTitleAbst.csv", row.names = F)
pubtatDatWithtTitlAbst <- pubtatDat[!grepl("[|]t|[|]a", pubtatDat)]
pubtatDatWithtTitlAbst <- pubtatDatWithtTitlAbst[pubtatDatWithtTitlAbst != ""]
df1<- as.data.frame(pubtatDatWithtTitlAbst)
pubtatDatWithtTitlAbst.df <- separate(df1, pubtatDatWithtTitlAbst, c("PMID","Start","End","Entity","Class","Entity_ID"), sep="\t")
write.csv(pubtatDatWithtTitlAbst.df, file = "pubtatEntity.csv", row.names = F)


# Create network for PMID--Entity with all the entity properties as node properties
net.df <- read.csv("pubtator/pubtator_central_export_pubtatEntity.csv", header = T, stringsAsFactors = F)
net.df <- net.df[,c(1,4,2,3,5,6)]
PMID_Ent_Net <- graph_from_data_frame(net.df, directed = FALSE, vertices = NULL)
# Write graph
write_graph(PMID_Ent_Net, "pubtator/PMID-Ent_Net.graphml", "graphml")

# Network creation between PMID and EntityID with class
net.df <- read.csv("pubtator/pubtator_central_export_pubtatEntity.csv", header = T, stringsAsFactors = F)
netMerz.df <- net.df[,c(1,4:6)]
netMerz.df <- unite(netMerz.df, EntityProp, c(Class, Entity_ID), remove=T)

for (i in 1:dim(netMerz.df)[1]) {
  avalEnts <- netMerz.df$Entity[which(netMerz.df$EntityProp %in% netMerz.df$EntityProp[i])]
  netMerz.df$Entity_CanonicalName[i] <- avalEnts[which(nchar(avalEnts) == min(nchar(avalEnts)))][1]
}

netMerz.df2 <- netMerz.df[c(1,3,2,4)]
PMID_Ent_NetId <- graph_from_data_frame(netMerz.df2, directed = FALSE, vertices = NULL)
# Write graph
write_graph(PMID_Ent_NetId, "pubtator/PMID-Ent_NetId.graphml", "graphml")


# Network creation between EntityID with class based on PMID
## Idea is to create edge between unique entity IDs of different classes if found in the same article
# Node file creation
## Columns: Node, Freq, CaptrdEnts
node.df <- data.frame(table(netMerz.df$EntityProp))
colnames(node.df)[1] <- "Node"
for (i in 1:dim(node.df)[1]) {
  entCaptrd <- node.df$Entity[which(netMerz.df$EntityProp %in% node.df$Node[i])]
  entCaptrd <- unique(entCaptrd)
  nodes$CaptrdEnts[i] <- paste(entCaptrd, collapse = ";")
}
write.csv(nodes, file = "pubtator/nodes.csv", row.names = F)

# Edge file creation
## Columns: Node1, Node2, PMID
df2 <- netMerz.df[,-2]
df2 <- unique(df2)
unikPMIDs <- unique(df2$PMID)
edgeList <- list()
for (i in 1:length(unikPMIDs)) {
  nodes <- df2$EntityProp[which(df2$PMID %in% unikPMIDs[i])]
  edge.df <- expand.grid(nodes, nodes)
  edge.df$PMID <- unikPMIDs[i]
  edge.df <- edge.df[!(edge.df$Var1 == edge.df$Var2),]
  edgeList[[i]] <- edge.df
  
}
edges = do.call(rbind, edgeList)
colnames(edges)[1:2] <- c("Node1", "Node2")
write.csv(edges, file = "pubtator/edges.csv", row.names = F)



# Creation of a network for top 10 genes with top gene i.e. Cancer
##Read the edge file for bioconcept-bioconcept network
setwd("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubtator/TopNetwrk")
edges <- read.csv("edges.csv", header = T, stringsAsFactors = F)
##Subset the file based on the nodes
nodes4Subset <- c("Disease_MESH:D009369","Gene_7157","Gene_7431","Gene_207","Gene_3569",
                  "Gene_7124","Gene_1956","Gene_335","Gene_2335","Gene_3315","Gene_2475")
edgesSubset <- edges[(edges$Node1 %in% nodes4Subset & edges$Node2 %in% nodes4Subset),]
##Change the node labels from ID into actual names
edgesSubset[edgesSubset == "Disease_MESH:D009369"] <- "Cancer"
edgesSubset[edgesSubset == "Gene_7157"] <- "p53"
edgesSubset[edgesSubset == "Gene_7431"] <- "Vimentin"
edgesSubset[edgesSubset == "Gene_207"] <- "AKT"
edgesSubset[edgesSubset == "Gene_3569"] <- "IL-6"
edgesSubset[edgesSubset == "Gene_7124"] <- "TNF-alpha"
edgesSubset[edgesSubset == "Gene_1956"] <- "EGFR"
edgesSubset[edgesSubset == "Gene_335"] <- "Apo A-I"
edgesSubset[edgesSubset == "Gene_2335"] <- "Fibronectin"
edgesSubset[edgesSubset == "Gene_3315"] <- "HSP27"
edgesSubset[edgesSubset == "Gene_2475"] <- "mTOR"
write.csv(edgesSubset, file = "edgesSubset.csv", row.names = F)
##Modify the edge file to contain weight and Gephi format
edgesSubsetModfd <- edgesSubset
###Add weight column
edgesSubsetModfd$Weight <- 1
edgesSubsetModfd <- aggregate(cbind(PMID, Weight) ~ Node1 + Node2, data=edgesSubsetModfd, function(x) cbind(paste(x, collapse = "|"), sum(as.numeric(x))))
edgesSubsetModfd$PMID <- edgesSubsetModfd$PMID[,1]
edgesSubsetModfd$Weight <- edgesSubsetModfd$Weight[,2]
edgesSubsetModfd$Types <- "Undirected"
colnames(edgesSubsetModfd)[1:2] <- c("Source", "Target")
edgesSubsetModfd <- edgesSubsetModfd[,c(1:2,5,4,3)]
write.csv(edgesSubsetModfd, file = "edgesSubset_Modfd4Gephi.csv", row.names = F)


# Create a Venn diagram or alternative for intersection between PMIDS under subheadings
##Alternative to Venn diagram; use upsetr; https://github.com/hms-dbmi/UpSetR/; https://gehlenborglab.shinyapps.io/upsetr/
##For this a file need to be prepared with 1st column as PMIDs and all subheadings as column names with 1 if PMID is and 0 if not

###merging the PMIDs to find unique all
setwd("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubmed/subheadings")
filelist <- list.files(pattern = ".*.txt")
datalist <- lapply(filelist, function(x)read.table(x, header=F))
allPMIDWthSubhdng.df <- do.call("rbind", datalist)
uniqPMIDWthSubhdng.df <- unique(allPMIDWthSubhdng.df)
colnames(uniqPMIDWthSubhdng.df) <- "PMID"
#uniqPMIDsUndSubhdngs <- unique(allPMID.df$V1)
###read and define all sets
agonists <- readLines("agonists_5.txt")
PMIDtemp <- c("12197311","31174591","22967081","31887191","20126400","29587743","30223649","26743025")
PMIDtemp.df <- as.data.frame(PMIDtemp)
colnames(PMIDtemp.df) <- "PMID"
PMIDtemp.df$agonists <- 0
for (i in 1:dim(PMIDtemp.df)[1]) {
  if(PMIDtemp.df$PMID[i] %in% agonists){
    PMIDtemp.df$agonists[i] <- 1
  }
}

###A function which reads all the files which are to be considered as subheadings and create a file with
###1st column as PMIDs and all subheadings as column names with 1 if PMID is and 0 if not
inpt4upsetr <- function(dir4subheadings){
  listFiles <- list.files(path=dir4subheadings, pattern = ".*.txt")
  datalist <- lapply(listFiles, function(x)read.table(x, header=F))
  allPMIDWthSubhdng.df <- do.call("rbind", datalist)
  uniqPMIDWthSubhdng.df <- unique(allPMIDWthSubhdng.df)
  colnames(uniqPMIDWthSubhdng.df) <- "PMID"
  ###read and define all sets
  for (i in 1:length(listFiles)) {
    fileNameTemp <- gsub("_(\\d+).txt","",listFiles[i])
    fileName <- gsub(" ", "_", fileNameTemp)
    PMIDsubhdng <- readLines(listFiles[i])
    uniqPMIDWthSubhdng.df$newCol <- 0
    for (j in 1:dim(uniqPMIDWthSubhdng.df)[1]) {
      if(uniqPMIDWthSubhdng.df$PMID[j] %in% PMIDsubhdng){
        uniqPMIDWthSubhdng.df$newCol[j] <- 1
      }
    }
    colnames(uniqPMIDWthSubhdng.df)[ncol(uniqPMIDWthSubhdng.df)] <- paste0(fileName)
  }
  write.csv(uniqPMIDWthSubhdng.df, file = paste0(dir4subheadings,"/","inptPMIDs4upsetr.csv"), row.names = F)
}

###Run
inpt4upsetr("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubmed/subheadings")



###
library(UpSetR)
inptPMIDs4upsetr <- read.csv("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubmed/temp/inptPMIDs4upsetr.csv", header = T)
upset(inptPMIDs4upsetr, sets = colnames(inptPMIDs4upsetr)[-1], sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

inptPMIDs4upsetr <- read.csv("inptPMIDs4upsetr.csv", header = T)
upset(inptPMIDs4upsetr, sets = colnames(inptPMIDs4upsetr)[-1], sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

# Modification of the edge file before importing onto Gephi
##There are some of the classes with no entity IDs, they can be removed
## CellLine_, Chemical_, Chemical_-, Gene_, Disease_

##Read the edge file
edgeFile <- read.csv("E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubtator/edges.csv", header = T, stringsAsFactors = F)
edgeFileModfd <- edgeFile[!(edgeFile$Node1 %in% c("CellLine_", "Chemical_", "Chemical_-", "Gene_", "Disease_")|edgeFile$Node2 %in% c("CellLine_", "Chemical_", "Chemical_-", "Gene_", "Disease_")),]
##Add weight column
edgeFileModfd$Weight <- 1
edgeFileModfd <- aggregate(cbind(PMID, Weight) ~ Node1 + Node2, data=edgeFileModfd, function(x) cbind(paste(x, collapse = "|"), sum(as.numeric(x))))
edgeFileModfd$PMID <- edgeFileModfd$PMID[,1]
edgeFileModfd$Weight <- edgeFileModfd$Weight[,2]
edgeFileModfd$Types <- "Undirected"
colnames(edgeFileModfd)[1:2] <- c("Source", "Target")
edgeFileModfd <- edgeFileModfd[,c(1:2,5,4,3)]
write.csv(edgeFileModfd, file = "E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubtator/edges_Modfd4Gephi.csv", row.names = F)

# SUbset for a particular bioconcept type
##For Gene--Gene associations
##using the above dataframe edgeFileModfd or else that file (or any file of interest) can be read and follow the same
###variable for bioconcept type
typeBioconcept1 <- "Gene"
typeBioconcept2 <- "Gene"
##create a match grep option to subset
mtchs <- grepl(paste0(typeBioconcept1,"_"), edgeFileModfd$Source, fixed=TRUE) & grepl(paste0(typeBioconcept2,"_"), edgeFileModfd$Target, fixed=TRUE)
edgeFileModfd4Gene <- edgeFileModfd[mtchs,]
write.csv(edgeFileModfd4Gene, file = "E:/Mine/ProteomesInvitatn/ProteomeArticleMinin/pubtator/edgesGene_Modfd4Gephi.csv", row.names = F)



###############################END###############################

