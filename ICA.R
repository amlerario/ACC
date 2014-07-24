#source("http://bioconductor.org/biocLite.R")
#biocLite(pkgs=c("Rsubread","limma"))
library(Biobase)
library(plyr)
library(ggplot2)
library(foreach)
library(xtable)
library(biomaRt)
library(GOstats)
library(cluster)
library(marray)
library(mclust)
library(RColorBrewer)
library(igraph)
library(Rgraphviz)
library(graph)
library(colorspace)
library(annotate)
library(scales)
library(gtools)
library(MineICA)
library(genefilter)
library(org.Hs.eg.db)
library(impute)
library(xlsx)

setwd("/home/amlerario/Projects/TCGA/EXPRESSION")


### edgeR.log2.txt
xx <- as.matrix(read.table("edgeR.log2.txt", header=TRUE, sep="\t", row.names=1, as.is=TRUE))
pData <- read.table("pData.txt", header=T,sep="\t", row.names=1, as.is=TRUE)
### REESCALING TO EXCLUDE NEGATIVE VALUES
#xxantilog <- apply(xx,2,function(x){ 2^x})
#constant<-1.001-range(as.vector(xxantilog), na.rm=T)[1]
#xxantilog <- apply(xxantilog,2,function(x){constant+x})
#xx <- apply(xxantilog,2,log2)
#xx[is.na(xx)] <- 1.001 # remove NA - attribute the minimal value to NA

## ANNOTATING
symbol<-unlist(lookUp(rownames(xx), 'org.Hs.eg', 'SYMBOL'))
symbol<-unname(symbol, force=FALSE)
name<-unlist(lookUp(rownames(xx), 'org.Hs.eg', 'GENENAME'))
name<-unname(name, force=FALSE)
assayData<-cbind(symbol,name)
rownames(assayData)<-symbol

assayData <- assayData[!is.na(rownames(assayData)),]
assayData <- data.frame(assayData)


rownames(xx)<-symbol
rownames(pData) <- colnames(xx)
xx <- xx[!is.na(rownames(xx)),]##remove those rows without SYMBOL

rownames(xx)<-as.factor(rownames(xx))
### - imputation - elliminate NAs from the table
#xx<-impute.knn(xx)

### create an eset object
phenoData <- new("AnnotatedDataFrame", data=pData)
assayData<-new("AnnotatedDataFrame", data=assayData)
#eset <- ExpressionSet(assayData=xx$data, phenoData=phenoData, annotation="SYMBOL",featureData=assayData)
eset <- ExpressionSet(assayData=xx, phenoData=phenoData, annotation="SYMBOL",featureData=assayData)
eset_var<-eset

library(genefilter)
exprs(eset_var) <- 2^exprs(eset_var)
#Filter #1 - exclude genes with very low expression values - consider 2.5 to 5 inferior percentile as "zero"
exprs(eset_var)[is.na(exprs(eset_var))] <- quantile(exprs(eset_var),0.01,na.rm=T) ## eliminate NA by substituting them by the p1 value
#eset_var <- selectFeatures_IQR(eset_var,10000)

ffun <- filterfun(pOverA(p=.8,A=quantile(exprs(eset_var),0.025,na.rm=T)))
#ffun <- filterfun(pOverA(p=0, A=-8))
t.fil <- genefilter(exprs(eset_var),ffun)
##################t.fil <- apply(e.mat,2,allNA)
# apply filter, and put expression back on log scale
exprs(eset_var) <- log2(exprs(eset_var)[t.fil,])
#zz<-impute.knn(exprs(eset_var))
#exprs(eset_var)<-zz$data

#ACC_RNASEQ_PRADA_RefGene
#xx <- as.matrix(read.table("ACC_RNASEQ_PRADA_RefGene.txt", header=TRUE, sep="\t", row.names=1, as.is=TRUE))
#isexpr <- rowSums(xx >= 0.5) >= length(xx[1,])*0.5
#xx <- xx[isexpr,]
#constant<-1.001-range(as.vector(xx), na.rm=T)[1]
#xx <- apply(xx,2,function(x){constant+x})
#xx[is.na(xx)] <- 1.001 # remove NA - attribute the minimal value to NA
#ffun <- filterfun(cv(0.7,10))
#t.fil <- genefilter(e.mat,ffun)
# apply filter, and put expression back on log scale
#xx <- log2(e.mat[t.fil,])
#xx <- apply(xx,2,log2)

## ANNOTATING PRADA
#xx <- xx[!is.na(rownames(xx)),]##remove those rows without SYMBOL
#rownames(pData) <- colnames(xx)


## Adding phenotype data and building an eSet object
#phenoData <- new("AnnotatedDataFrame", data=pData)
#eset <- ExpressionSet(assayData=xx, phenoData=phenoData, annotation="hgnc_symbol")


# obtain expression estimates on the UN-LOGGED scale

#e.mat <- 2^exprs(eset)
# look at mean, sd, & cv for each gene across arrays
#gene.mean <- apply(e.mat,1,mean)
#gene.sd <- apply(e.mat,1,sd)
#gene.cv <- gene.sd/gene.mean
# make plots
#library(geneplotter); library(RColorBrewer)
#blues.ramp <- colorRampPalette(brewer.pal(9,"Blues")[-1])
#dCol <- densCols(log(gene.mean),log(gene.sd),colramp=blues.ramp)
#par(mfrow=c(2,2))
#plot(gene.mean,gene.sd,log='xy',col=dCol,pch=16,cex=0.1)
#abline(v=5,lwd=3,col='red')
#hist(log(gene.cv),main=NA)
#abline(v=log(.7),lwd=3,col='red')

## Filtering by variance
#eset_var <- selectFeatures_IQR(eset,round(dim(exprs(eset))[1]*0.5)) # IQR p50

## restrict the number of genes to 10000 (as suggested by Anne Biton)
#eset_var <- selectFeatures_IQR(eset_var,10000)

########################################################################
## D E F I N E  H E R E  T H E  C U T O F F  V A L U E  F O R  I Q R  ##
cutoff_iqr=0   #########################################################
########################################################################

#if(cutoff_iqr != 0) {
#  eset_var <- nsFilter(eset_var, var.cutoff=cutoff_iqr/100, require.entrez=F, remove.dupEntrez=F)
#  eset_var <- eset_var$eset
#}


library(JADE)

## Features are mean-centered before ICA computation

exprs(eset_var) <- t(apply(exprs(eset_var),1,scale,scale=FALSE))#scale=FALSE
colnames(exprs(eset_var)) <- sampleNames(eset_var)

## Run ICA-JADE

########################################################################
## D E F I N E  H E R E  T H E  N U M B E R  O F  C O M P O N E N T S ##
ncomp=8  ###############################################################
########################################################################

########################################################################
##              D E F I N E  H E R E  T H E  C U T O F F S            ##
selCutoff=3.5  #########################################################
########################################################################

resJade <- runICA(X=exprs(eset_var), nbComp=ncomp, method="JADE", maxit=dim(exprs(eset_var))[1])
## if an error message appears here, try to reduce the ncomp variable
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
typeIDeset_var <-c(geneID_annotation="SYMBOL", geneID_biomart="hgnc_symbol")#"hgnc_symbol"entrezgene
params <- buildMineICAParams(resPath=paste("ICA.com",ncomp,"_iqr",cutoff_iqr,"_",selCutoff, "/",sep=""), selCutoff=selCutoff, pvalCutoff=0.01)
refSamplesMainz <- character(0)
resBuild <- buildIcaSet(params=params, A=data.frame(resJade$A), S=data.frame(resJade$S), dat=exprs(eset_var),
                        pData=pData(eset_var),refSamples=refSamplesMainz, typeID=typeIDeset_var,mart=mart)
icaSeteset_var <- resBuild$icaSet
params <- resBuild$params

## see sample projections of specific components - in this example, component 1
#comp1<-getComp(icaSeteset_var, level="genes", ind=1)
#comp18<-getComp(icaSeteset_var, level="genes", ind=18)

keepVar <- c("TP53","CTNNB1","ZNRF3_ALT","weiss","Mitosis_50","Grade","hormone_excess","histological_type","age_diagnosis","ENSAT", "mRNA_K5", "mRNA_K2", "Methy_CLS", "miRNA_CLS", "SCNA_CLS", "LOH", "n_mutation", "gender", "vital_status", "PKA")
runAn(params=params, icaSet=icaSeteset_var, writeGenesByComp=TRUE, keepVar=keepVar,mart=mart,selCutoffWrite=selCutoff)

#resW <- writeProjByComp(icaSet=icaSeteset_var,params=params,mart=mart,level="genes",selCutoffWrite=3)

#[1] "mRNA_K5"               "mRNA_K2"               "ENSAT"                 "Methy_CLS"            
#[5] "MethyLevel_CLS"        "miRNA_CLS"             "SCNA_CLS"              "PKA"                  
#[9] "LOH"                   "purity"                "ploidy"                "n_mutation"           
#[13] "hypermut"              "histological_type"     "gender"                "vital_status"         
#[17] "weiss"                 "Mitosis_50"            "hormone_excess"        "age_diagnosis"        
#[21] "days_to_birth"         "days_to_death"         "days_to_last_followup" "TP53_mut"             
#[25] "TP53_SCNA"             "CTNNB1_mut"            "CTNNB1_SCNA"           "MEN1_mut"             
#[29] "MEN1_SCNA"             "PRKAR1A_mut"           "PRKAR1A_SCNA"          "RPL22_mut"            
#[33] "RPL22_SCNA"            "X12q14_CDK4"           "X16q22_TERF2"          "X17q25_RFNG"          
#[37] "X19p13_BRD4_.AKAP8"    "X19q12_CCNE1"          "X1q22_EFNA3_4"         "X5p15_TERT"           
#[41] "xq28_PNMA6A"           "X22q12_ZNRF3"          "ZNRF3_mut"             "X9p21_CDKN2A"         
#[45] "CDKN2A_methelation"    "APC"                   "ZNRF3"                 "TP53" 


library(STRINGdb)
string_db <- STRINGdb$new(version="9_1", species=9606, score_threshold=0, input_directory="")
contrib <- selectContrib(icaSeteset_var, cutoff=3, level="genes")
#genes18<-sort(contrib[[18]],decreasing=TRUE)#[1:10]

#genes18<-data.frame(names(genes18))

### evaluating interactions of the positive side of the components

dir.create(paste(getwd(),"/ICA.com",ncomp,"_iqr",cutoff_iqr,"_",selCutoff,"/string",sep=""),showWarnings=F)

for(i in 1:ncomp) {
  genes<-sort(contrib[[i]][contrib[[i]]>0],decreasing=TRUE)
  genes<-data.frame(names(genes))
  mapped <- string_db$map(genes,"names.genes.", removeUnmappedRows = TRUE)
  hits <- mapped$STRING_id
  string_db$plot_network(hits)
  dev.copy(pdf,paste(getwd(),"/ICA.com",ncomp,"_iqr",cutoff_iqr,"_",selCutoff,"/string/","comppos",i,".pdf",sep=""))
  dev.off()
}

### evaluating interactions of the negative side of the components
### create 

for(i in 1:ncomp) {
  genes<-sort(contrib[[i]][contrib[[i]]<0],decreasing=TRUE)
  genes<-data.frame(names(genes))
  mapped <- string_db$map(genes,"names.genes.", removeUnmappedRows = TRUE)
  hits <- mapped$STRING_id
  string_db$plot_network(hits)
  dev.copy(pdf,paste(getwd(),"/ICA.com",ncomp,"_iqr",cutoff_iqr,"_",selCutoff,"/string/","compneg",i,".pdf",sep=""))
  dev.off()
}


for(i in 1:ncomp) {
  genesneg<-sort(contrib[[i]][contrib[[i]]<0],decreasing=TRUE)
  genesneg<-data.frame(names(genesneg))
  mappedneg <- string_db$map(genesneg,"names.genesneg.", removeUnmappedRows = TRUE)
  hitsneg <- mappedneg$STRING_id
  enrichmentGOneg <- string_db$get_enrichment(hitsneg, category = "Process", methodMT = "fdr", iea = TRUE)
  enrichmentKEGGneg <- string_db$get_enrichment(hitsneg, category = "KEGG", methodMT = "fdr", iea = TRUE )
  genespos<-sort(contrib[[i]][contrib[[i]]>0],decreasing=TRUE)
  genespos<-data.frame(names(genespos))
  mappedpos <- string_db$map(genespos,"names.genespos.", removeUnmappedRows = TRUE)
  hitspos <- mappedpos$STRING_id
  enrichmentGOpos <- string_db$get_enrichment(hitspos, category = "Process", methodMT = "fdr", iea = TRUE)
  enrichmentKEGGpos <- string_db$get_enrichment(hitspos, category = "KEGG", methodMT = "fdr", iea = TRUE )
  genes<-sort(contrib[[i]],decreasing=TRUE)
  genes<-data.frame(names(genes))
  mapped <- string_db$map(genes,"names.genes.", removeUnmappedRows = TRUE)
  hits <- mapped$STRING_id
  enrichmentGO <- string_db$get_enrichment(hits, category = "Process", methodMT = "fdr", iea = TRUE)
  enrichmentKEGG <- string_db$get_enrichment(hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
  for(j in 1:6){
    if(j==1){
      sheet<-"GO_NEG"
      data <- enrichmentGOneg 
    }
    if(j==2){
      sheet<-"KEGG_NEG"
      data <- enrichmentKEGGneg 
    }
    if(j==3){
      sheet<-"GO_POS"
      data <- enrichmentGOpos 
    }
    if(j==4){
      sheet<-"KEGG_POS"
      data <- enrichmentKEGGpos 
    }
    if(j==5){
      sheet<-"GO_ALL"
      data <- enrichmentGO 
    }
    if(j==6){
      sheet<-"KEGG_ALL"
      data <- enrichmentKEGG 
    }
    write.xlsx(data.frame(data),file=paste(getwd(),"/ICA.com",ncomp,"_iqr",cutoff_iqr,"_",selCutoff,"/string/","comp",i,".xlsx",sep=""),sheetName=sheet,col.names=T,row.names=T,append=T,showNA=T)
  }
  
}



#geneint=c("GATA4")
#closeg <- genefinder(exprs(eset), geneint, 200, method="maximum", scale="zscore")
#rownames(exprs(eset))[closeg[[1]]$indices]


#resEnrich <- runEnrich(params=params,icaSet=icaSeteset_var[,,3],dbs="GO", ontos="BP", cond=T)

#head(icaSeteset_var[,1:20])

### WTINESS GENES OF EACH COMPONENT

#witGenes(icaSeteset_var)
#1         2         3         4         5         6         7         8         9        10        11        12        13 
#"RPS4Y1"   "GSTA1"  "SPRR1A"     "NDN"  "CTNNA2"  "DLGAP5"  "NKAIN4"     "LUM" "TDGF1P3" "RASGRF1"   "PROK1"    "MAOB" "FAM166B" 
#14        15        16        17        18        19        20        21        22 
#"B3GNT7"   "FCGBP" "FAM19A5"    "ZFR2" "SLC4A10"    "ORM1"   "EGLN3"   "MGAT5"    "SNCB" 