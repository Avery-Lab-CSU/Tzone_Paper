library(AnnotationDbi)
library(org.Cf.eg.db)
library(biomaRt)
library(msigdbr)
library(enrichplot)
library(vroom)
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(magrittr)
library(data.table)
# Make raw count matrix
setwd("~/Documents/Avery_lab/tzn_100annotation/remove1/cnts//")
filelist=Sys.glob("C*")
f1_counts<-read.table(filelist[1],sep='\t')
number_rows<-nrow(f1_counts)
raw_counts<-data.frame(f1_counts[5:number_rows,2])
rownames(raw_counts)<-f1_counts[5:number_rows,1]
colnames(raw_counts)[1]<-substr(filelist[1],start = 1, stop = nchar(filelist[1])-21)	
number_files<-length(filelist)
for (file in filelist[2:number_files]){
  holder<-read.table(file,row.names=1, sep='\t')
  colnames(holder)[1]<-substr(file,start = 1, stop = nchar(file)-21)	
  raw_counts<-merge(raw_counts,holder,by='row.names')
  rownames(raw_counts)<-raw_counts[,1]
  raw_counts<-raw_counts[,-1]
  number_columns<-ncol(raw_counts)
  raw_counts<-raw_counts[,-(number_columns)]
  raw_counts<-raw_counts[,-(number_columns-1)]
}
filelist=Sys.glob("T*")
for (file in filelist){
  holder<-read.table(file,row.names=1, sep='\t')
  colnames(holder)[3]<-substr(file,start = 1, stop = nchar(file)-21)	
  raw_counts<-merge(raw_counts,holder,by='row.names')
  rownames(raw_counts)<-raw_counts[,1]
  raw_counts<-raw_counts[,-1]
  number_columns<-ncol(raw_counts)
  raw_counts<-raw_counts[,-c((number_columns-2),(number_columns-1))]
}
raw_counts=raw_counts[,c(2:4,1,5:13)]

# Metadata
cd8cd4md=data.frame(condition=c(rep.int("cd8_cntrl",3),rep.int("cd4_control",4)), 
                    row.names=colnames(raw_counts[,1:7]))
tzncd8md=data.frame(condition=c(rep.int("cd8_cntrl",3),rep.int("tzn",6)), 
                    row.names=colnames(raw_counts[,c(1:3,8:13)]))
cd4md=data.frame(condition=c(rep.int("cd4_cntrl",4),rep.int("tzn",6)), 
                 row.names=colnames(raw_counts[,c(4:7,8:13)]))
alldmd=data.frame(condition=c(rep.int("cd8_cntrl",3),
                              rep.int("cd4_cntrl",4),
                              rep.int("tzn",6)), 
                  row.names=colnames(raw_counts))

#DeSeq data sets and run
cd8cd4dds=DESeqDataSetFromMatrix(countData = raw_counts[,1:7], 
                                 colData = cd8cd4md, 
                                 design = ~ condition)
cd8cd4dds<-DESeq(cd8cd4dds)
cd4dds=DESeqDataSetFromMatrix(countData = raw_counts[,c(4:13)], 
                              colData = cd4md, 
                              design = ~ condition)
tzncg8dds=DESeqDataSetFromMatrix(countData = raw_counts[,c(1:3,8:13)], 
                                 colData = tzncd8md, 
                                 design = ~ condition)
cd4dds<-DESeq(cd4dds)
tzncg8dds=DESeq(tzncg8dds)
alldds=DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = alldmd, 
                              design = ~ condition)
alldds=DESeq(alldds)

#Get result dataframes
res3=results(tzncg8dds, tidy=T)
res=results(cd4dds, tidy = T)
cres=results(cd8cd4dds, tidy = T)
rownames(res)=res[,1]
rownames(res3)=res3[,1]
rownames(cres)=cres[,1]
res=res[,-1]
res3=res3[,-1]
cres=cres[,-1]

#Get normalized counts, rearrange columns
normc=counts(alldds, normalized=T)
finalres=merge(normc, res3, by="row.names")
rownames(finalres)=finalres[,1]
finalres=finalres[,-1]
finalres=merge(finalres, res, by="row.names")
rownames(finalres)=finalres[,1]
finalres=finalres[,-1]
finalres=finalres[,c(1:13,19,25)]
colnames(finalres)[14]="FDR_tzn_cd8"
colnames(finalres)[15]="FDR_tzn_cd4"

#Calc means
finalres$cd8_mean=rowMeans(finalres[,1:3])
finalres$cd4_mean=rowMeans(finalres[,4:7])
finalres$tzn_mean=rowMeans(finalres[,8:13])

#If mean is zero make 1 for fold change
finalres[finalres$cd8_mean==0,16]=1
finalres[finalres$cd4_mean==0,17]=1
finalres[finalres$tzn_mean==0,18]=1

finalres$l2fc_tzn_cd8=log2(finalres$tzn_mean / finalres$cd8_mean)
finalres$l2fc_tzn_cd4=log2(finalres$tzn_mean / finalres$cd4_mean)
finalres$l2fc_cd8_cd4=log2(finalres$cd8_mean / finalres$cd4_mean)

finalres=merge(finalres, cres, by="row.names")
rownames(finalres)=finalres[,1]
finalres=finalres[,-1]
colnames(finalres)[27] <- "FDR_cd8_cd4"

finalres=finalres[,c(1:15,27,16:21)]

# Configure biomaRt for archive version 100
dog = useMart("ensembl", host = "https://apr2020.archive.ensembl.org",
              dataset = "clfamiliaris_gene_ensembl")
dogbm=getBM(attributes = c("ensembl_gene_id", 
                           "hgnc_symbol", 
                           "uniprot_gn_symbol",
                           "entrezgene_id"), mart = dog)

dogbm2=getBM(attributes = c("ensembl_gene_id", 
                            "description",
                            "entrezgene_id",
                            "entrezgene_accession"), mart = dog)

# Functions for converting eids
dens2hgncsym=function(ens){
  idx1 = which(dogbm$ensembl_gene_id == ens)
  if(length(idx1) == 0){return(NA)
  }else{
    sym = dogbm$hgnc_symbol[idx1[1]]
    return(sym)
  }
}

dens2unisym=function(ens){
  idx1 = which(dogbm$ensembl_gene_id == ens)
  if(length(idx1) == 0){return(NA)
  }else{
    sym = dogbm$uniprot_gn_symbol[idx1[1]]
    return(sym)
  }
}

dens2entid=function(ens){
  idx1 = which(dogbm2$ensembl_gene_id == ens)
  if(length(idx1) == 0){return(NA)
  }else{
    ent = dogbm2$entrezgene_accession[idx1[1]]
    return(ent)
  }
}

dens2entacc=function(ens){
  idx1 = which(dogbm2$ensembl_gene_id == ens)
  if(length(idx1) == 0){return(NA)
  }else{
    ent = dogbm2$entrezgene_id[idx1[1]]
    return(ent)
  }
}

dens2desc=function(ens){
  idx1 = which(dogbm2$ensembl_gene_id == ens)
  if(length(idx1) == 0){return(NA)
  }else{
    ent = dogbm2$description[idx1[1]]
    return(ent)
  }
}

# Convert 
finalres$uniprot_sym=sapply(rownames(finalres), FUN = dens2unisym, simplify = T, USE.NAMES = F)
finalres$hgnc_sym=sapply(rownames(finalres), FUN = dens2hgncsym, simplify = T, USE.NAMES = F)
finalres$entrez <- sapply(rownames(finalres), FUN = dens2entid, simplify = T, USE.NAMES = F)
finalres$entrezid <- sapply(rownames(finalres), FUN = dens2entacc, simplify = T, USE.NAMES = F)
finalres$description <- sapply(rownames(finalres), FUN = dens2desc, simplify = T, USE.NAMES = F)

# Rename columns
for (i in 1:3){
  colnames(finalres)[i]<-paste("CD8",colnames(finalres)[i], sep="_")
}
for (i in 4:7){
  colnames(finalres)[i]<-paste("CD4",colnames(finalres)[i], sep="_")
}
cnames <- colnames(finalres)
finalres <- rownames_to_column(finalres)
colnames(finalres) <- c("ensembl_id",cnames)
rm(cnames, i, file, dogbm2, cd8cd4dds, cd8cd4md, cres, alldds, alldmd, cd4dds, cd4md, dog, dogbm, f1_counts,holder, normc, raw_counts, res, res3,tzncd8md, tzncg8dds, number_columns, number_files, number_rows, filelist)


########GSEA##########

sigres=finalres[finalres$FDR_tzn_cd8<.05,]
sigres=sigres[!is.na(sigres$FDR_tzn_cd8),]
sigres <- sigres[abs(sigres$l2fc_tzn_cd8) > 1.5,]
sigres <- sigres[!is.na(sigres$entrez),]

s2nReturn=function (expmat, label) {
  freq <- table(label)
  if ((freq[1] < 2) || (freq[2] < 2)) {
    stop("There are not enough samples for calculating singal to noise ratio ya dingdong!")
  }
  x0 <- expmat[, which(label == 0)]
  x1 <- expmat[, which(label == 1)]
  m0 <- apply(x0, 1, mean)
  m1 <- apply(x1, 1, mean)
  sd0 <- apply(x0, 1, sd)
  sd1 <- apply(x1, 1, sd)
  s2n <- (m1 - m0)/(sd0 + sd1)
  return(s2n)
}
sigres$s2n=s2nReturn(sigres[,c(2:4,9:14)],c(0,0,0,1,1,1,1,1,1))
geneList=sigres$s2n
names(geneList)=as.character(sigres$entrezid)
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!duplicated(names(geneList))]
m_t2g <- msigdbr(species = "Canis lupus familiaris", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
em <- GSEA(geneList, 
           TERM2GENE = m_t2g, 
           eps=0)

em2=setReadable(em, 'org.Cf.eg.db', 'ENTREZID')
em3 = filter(em2, !(Description %like% "KO"))
em3 = filter(em3, !(Description %like% "_DC_"))
em3 = filter(em3, !(Description %like% "_TRANSDUCED_"))
em3 = filter(em3, !(Description %like% "_BCELL_"))
em3 = filter(em3, !(Description %like% "_MONOCYTE_"))
em3 = filter(em3, !(Description %like% "_NKTCELL_"))
em3 = filter(em3, !(Description %like% "_DIABETIC_"))
em3 = filter(em3, !(Description %like% "_NKT_"))
em3 = filter(em3, !(Description %like% "DELETED"))
em3 = filter(em3, !(Description %like% "SPLENOCYTES"))
em3 = filter(em3, !(Description %like% "BMDC"))
em3 = filter(em3, !(Description %like% "MACROPHAGE"))
em3 = filter(em3, !(Description %like% "_DEFICIENT")) 
em3 = filter(em3, !(Description %like% "MONOCYTE"))
em3 = filter(em3, !(Description %like% "KNOCKIN_"))
em3 = filter(em3, !(Description %like% "NIH3T3"))
em3 = filter(em3, !(Description %like% "FUSION"))
em3 = filter(em3, !(Description %like% "MELANOMA"))
em3 = filter(em3, !(Description %like% "EPITHELIAL"))
em3 = filter(em3, !(Description %like% "_EOSINOPHIL_"))
em3 = filter(em3, !(Description %like% "ERYTHROD"))
em3 = filter(em3, !(Description %like% "MYELOD"))
em3 = filter(em3, !(Description %like% "_MICROGLIA_"))
em3 = filter(em3, !(Description %like% "MAST"))
em3 = filter(em3, !(Description %like% "IMPLANT"))
em3 = filter(em3, !(Description %like% "MAC"))
em3 = filter(em3, !(Description %like% "_PDC_"))
em3 = filter(em3, !(Description %like% "_EOSINOPHIL_"))
em3 = filter(em3, !(Description %like% "_NKCELL"))
em3 = filter(em3, !(Description %like% "_BMP_"))
em3 = filter(em3, !(Description %like% "MEF"))
em3 = filter(em3, !(Description %like% "_NEUTROPHIL"))
em3 = filter(em3, !(Description %like% "PBMC"))
em3 = filter(em3, !(Description %like% "BONE_MARROW"))
em3 = filter(em3, !(Description %like% "IRES"))

em3 %>%
  dotplot(showCategory = 49, x = "NES") +
  scale_colour_viridis_c(name = "Adjusted\nP-value",
                         option="B") +
  geom_vline(xintercept = 0, linetype=2)
dev.off()
