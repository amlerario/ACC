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
pData <- read.table("/home/amlerario/Projects/TCGA/EXPRESSION/scripts/pData.txt", header=T,sep="\t", row.names=1, as.is=TRUE)
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


## Run ICA-JADE



fish <- function(x,y) {
  table<-table(x,y)
  test<-fisher.test(table,workspace=2e8)
  print(table)
  print(test)
}

fish(categorical_variable1, categorical variable2)

### Survival analyses
## Kaplan-Meier

library(survival)
#attach(pData)

### evaluating factors that predict OS
for(tit in c("mRNA_K2","mRNA_K5","ENSAT","SCNA_CLS", "Methy_CLS", "MethyLevel_CLS","miRNA_CLS","TP53", "CTNNB1", "ZNRF3_ALT",
             "STAGE_AJCC","Grade", "hormone_excess", "cortisol_status")){
  title<-paste(tit,"_OS",sep="")
  sub_set <- pData
  fact <- sub_set[,tit]
  survdiff(Surv(sub_set[,"final_time"], sub_set[,"obs_death"])~fact, data=sub_set)
  result<-capture.output(survdiff(Surv(sub_set[,"final_time"], sub_set[,"obs_death"])~fact, data=sub_set))
  cat(c("########",title,result),file=paste(getwd(),"/survival/results.txt",sep=""),sep="\n",append=T)
  plot(survfit(Surv(sub_set[,"final_time"], sub_set[,"obs_death"])~fact, data=sub_set), main=title, lty=c(1,2,3,4,5,6), ylab="Prob.", 
       xlab="Survival time (days)", xlim=c(0,max(sub_set[,"final_time"], na.rm=T)*1.05))
  legend(200, 0.5, legend=levels(factor(fact)),lty=c(1,2,3,4,5,6), title=title, bty="n", cex=0.7)
  dev.copy(pdf,paste(getwd(),"/survival/",title,".pdf",sep=""))
  dev.off() 
  #### DFS
  title<-paste(tit,"_DFS",sep="")
  sub_set <- subset(pData, STAGE_AJCC != "Stage IV")
  fact<-sub_set[,tit]
  survdiff(Surv(sub_set[,"recurrence_time"], sub_set[,"recurrence"])~fact, data=sub_set)
  result<-capture.output(survdiff(Surv(sub_set[,"recurrence_time"], sub_set[,"recurrence"])~fact, data=sub_set))
  cat(c("########",title,result),file=paste(getwd(),"/survival/results_DFS.txt",sep=""),sep="\n",append=T)
  plot(survfit(Surv(sub_set[,"recurrence_time"], sub_set[,"recurrence"])~fact, data=sub_set), main=title, lty=c(1,2,3), ylab="Prob.", 
       xlab="Survival time (days)", xlim=c(0,max(sub_set[,"final_time"], na.rm=T)*1.05))
  legend(200, 0.5, legend=levels(factor(fact)),lty=c(1,2,3), title=title, bty="n", cex=0.7)
  dev.copy(pdf,paste(getwd(),"/survival/",title,"_DFS.pdf",sep=""))
  dev.off()
}




#### DFS
title<-paste("ZNRF3_ALT","_DFS", sep="")
sub_set <- subset(pData, STAGE_AJCC != "Stage IV")
fact<-sub_set$ZNRF3_ALT

survdiff(Surv(recurrence_time, recurrence)~fact, data=sub_set)
result<-capture.output(survdiff(Surv(recurrence_time, recurrence)~fact, data=sub_set))
cat(c("########",title,result),file=paste(getwd(),"/survival/results_DFS.txt",sep=""),sep="\n",append=T)
plot(survfit(Surv(recurrence_time, recurrence)~fact, data=sub_set), main=title, lty=c(1,2,3), ylab="Prob.", 
     xlab="Survival time (days)", xlim=c(0,max(final_time, na.rm=T)*1.05))
legend(200, 0.5, legend=levels(factor(fact)),lty=c(1,2,3), title=title, bty="n", cex=0.7)
dev.copy(pdf,paste(getwd(),"/survival/",title,"_DFS.pdf",sep=""))
dev.off()
