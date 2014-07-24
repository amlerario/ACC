source("http://bioconductor.org/biocLite.R")
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


#biocLite("Heatplus")
library(Heatplus)
#library(devtools)
#install_github("rasmusab/bayesian_first_aid")
library(BayesianFirstAid)

#install.packages("rjags")
#install.packages("mvtnorm")
#install.packages("car")
library(rjags)
library(mvtnorm)
library(car)


setwd("/home/amlerario/Projects/TCGA/EXPRESSION")

### edgeR.log2.txt
xx <- as.matrix(read.table("edgeR.log2.txt", header=TRUE, sep="\t", row.names=1, as.is=TRUE))
pData <- read.table("pData.txt", header=T,sep="\t", row.names=1, as.is=TRUE)
### REESCALING TO EXCLUDE NEGATIVE VALUES
xxantilog <- apply(xx,2,function(x){ 2^x})
constant<-1.001-range(as.vector(xxantilog), na.rm=T)[1]
xxantilog <- apply(xxantilog,2,function(x){constant+x})
xx <- apply(xxantilog,2,log2)

#xx[is.na(xx)] <- 1.001 # remove NA - attribute the minimal value to NA
## ANNOTATING
symbol<-unlist(lookUp(rownames(xx), 'org.Hs.eg', 'SYMBOL'))
symbol<-unname(symbol, force=FALSE)
rownames(xx)<-symbol
rownames(pData) <- colnames(xx)
xx <- xx[!is.na(rownames(xx)),]##remove those rows without SYMBOL

# filter: keep genes with cv between .7 and 10,
#
#and where 20% of samples had exprs. > 5
library(genefilter)
e.mat <- 2^xx
ffun <- filterfun(pOverA(0.2,5), cv(0.7,10))
t.fil <- genefilter(e.mat,ffun)
# apply filter, and put expression back on log scale
xx <- log2(e.mat[t.fil,])



#ACC_RNASEQ_PRADA_RefGene
xx <- as.matrix(read.table("ACC_RNASEQ_PRADA_RefGene.txt", header=TRUE, sep="\t", row.names=1, as.is=TRUE))
isexpr <- rowSums(xx >= 0.5) >= length(xx[1,])*0.5
xx <- xx[isexpr,]
constant<-1.001-range(as.vector(xx), na.rm=T)[1]
xx <- apply(xx,2,function(x){constant+x})
xx[is.na(xx)] <- 1.001 # remove NA - attribute the minimal value to NA
ffun <- filterfun(cv(0.7,10))
t.fil <- genefilter(e.mat,ffun)
# apply filter, and put expression back on log scale
xx <- log2(e.mat[t.fil,])
#xx <- apply(xx,2,log2)

## ANNOTATING
xx <- xx[!is.na(rownames(xx)),]##remove those rows without SYMBOL
rownames(pData) <- colnames(xx)


## Adding phenotype data and building an eSet object

xx<-impute.knn(xx)
phenoData <- new("AnnotatedDataFrame", data=pData)
eset <- ExpressionSet(assayData=xx$data, phenoData=phenoData, annotation="hgnc_symbol")
                      



# obtain expression estimates on the UN-LOGGED scale

#e.mat <- 2^exprs(eset)
# look at mean, sd, & cv for each gene across arrays
gene.mean <- apply(e.mat,1,mean)
gene.sd <- apply(e.mat,1,sd)
gene.cv <- gene.sd/gene.mean
# make plots
library(geneplotter); library(RColorBrewer)
blues.ramp <- colorRampPalette(brewer.pal(9,"Blues")[-1])
dCol <- densCols(log(gene.mean),log(gene.sd),
                 colramp=blues.ramp)
par(mfrow=c(2,2))
plot(gene.mean,gene.sd,log='xy',col=dCol,pch=16,cex=0.1)
abline(v=5,lwd=3,col='red')
hist(log(gene.cv),main=NA)
abline(v=log(.7),lwd=3,col='red')

## Filtering by variance
#eset_var <- selectFeatures_IQR(eset,round(dim(exprs(eset))[1]*0.5)) # IQR p50

eset_var <- eset


ttestfun <- function(z) {
  CTNNB1mut=z[eset$CTNNB1!="WT"]
  CTNNB1wt=z[eset$CTNNB1=="WT"]
  meanmut=mean(z[eset$CTNNB1!="WT"])
  meanwt=mean(z[eset$CTNNB1=="WT"])
  p=t.test(CTNNB1mut, CTNNB1wt, var.equal=T)$p.value
  foldchange=mean(CTNNB1mut)/mean(CTNNB1wt)
  c(meanmut,meanwt,foldchange,p)
}


results=t(apply(exprs(eset),1,ttestfun))
colnames(results) <- c("exprs_mut","exprs_wt","FC","p")
results <- results[order(results[,4]),]

BH=p.adjust(results[,4],method="fdr")
results=cbind(results,BH)

head(results,300)

library(gplots)
selected <- exprs(eset)[head(names(results[,1]),75),]
heatmap(selected, col=greenred(75), xlab="Conditions", ylab="Probes")


ttestfunWNT <- function(z) {
  WNTmut=z[eset$CTNNB1!="WT" | eset$ZNRF3!="WT"]
  WNTwt=z[eset$CTNNB1=="WT" &  eset$ZNRF3=="WT" ]
  meanmut=mean(z[eset$CTNNB1!="WT" | eset$ZNRF3!="WT"])
  meanwt=mean(z[eset$CTNNB1=="WT" &  eset$ZNRF3=="WT" ])
  p=t.test(WNTmut, WNTwt, var.equal=T)$p.value
  foldchange=mean(WNTmut)/mean(WNTwt)
  c(meanmut,meanwt,foldchange,p)
}

resultsWNT=t(apply(exprs(eset),1,ttestfunWNT))
colnames(resultsWNT) <- c("exprs_mut","exprs_wt","FC","p")
resultsWNT <- resultsWNT[order(resultsWNT[,4]),]

BH=p.adjust(resultsWNT[,4],method="fdr")
resultsWNT=cbind(resultsWNT,BH)

head(resultsWNT,300)

library(gplots)
selectedWNT <- exprs(eset)[head(names(resultsWNT[,1]),75),]
heatmap(log2(selectedWNT), col=greenred(75), xlab="Conditions", ylab="Probes")




# Writing tables
## by type
### contrast 1
write.table(head(results,300), file="PKAmut-PKAwt.xls", row.names=TRUE, col.names=TRUE, sep="\t")

boxplot(as.vector(exprs(eset[rownames(eset)=="LEF1",]))~eset$WNT, data=exprs(eset), main="Expression levels of LEF1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="LEF1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="AXIN2",]))~eset$WNT, data=exprs(eset), main="Expression levels of AXIN2 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="AXIN2",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="TCF7",]))~eset$WNT, data=exprs(eset), main="Expression levels of TCF7 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="TCF7",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="ZNRF3",]))~eset$WNT, data=exprs(eset), main="Expression levels of ZNRF3 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="ZNRF3",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="LGR5",]))~eset$WNT, data=exprs(eset), main="Expression levels of LGR5 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="LGR5",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="BMP4",]))~eset$WNT, data=exprs(eset), main="Expression levels of BMP4 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="BMP4",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="NOV",]))~eset$WNT, data=exprs(eset), main="Expression levels of NOV according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="NOV",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="CCDC80",]))~eset$WNT, data=exprs(eset), main="Expression levels of CCDC80 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="CCDC80",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="TNFSF4",]))~eset$WNT, data=exprs(eset), main="Expression levels of TNFSF4 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="TNFSF4",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="ISM1",]))~eset$WNT, data=exprs(eset), main="Expression levels of ISM1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="ISM1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="NKD1",]))~eset$WNT, data=exprs(eset), main="Expression levels of NKD1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="NKD1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="NR5A1",]))~eset$WNT, data=exprs(eset), main="Expression levels of NR5A1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="NR5A1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="CYP11B2",]))~eset$WNT, data=exprs(eset), main="Expression levels of CYP11B2 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="CYP11B2",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="KCNJ5",]))~eset$WNT, data=exprs(eset), main="Expression levels of KCNJ5 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="KCNJ5",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="IGF2",]))~eset$WNT, data=exprs(eset), main="Expression levels of IGF2 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="IGF2",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="CTNNB1",]))~eset$WNT, data=exprs(eset), main="Expression levels of CTNNB1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="CTNNB1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="RNF43",]))~eset$WNT, data=exprs(eset), main="Expression levels of RNF43 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="RNF43",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="CCND1",]))~eset$WNT, data=exprs(eset), main="Expression levels of CCND1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="CCND1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="MYC",]))~eset$WNT, data=exprs(eset), main="Expression levels of MYC according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="MYC",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="MYB",]))~eset$WNT, data=exprs(eset), main="Expression levels of MYB according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="MYB",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="TOP2A",]))~eset$WNT, data=exprs(eset), main="Expression levels of TOP2A according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="TOP2A",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="APC",]))~eset$WNT, data=exprs(eset), main="Expression levels of APC according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="APC",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="CSGALNACT1",]))~eset$WNT, data=exprs(eset), main="Expression levels of CSGALNACT1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="CSGALNACT1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="RELN",]))~eset$WNT, data=exprs(eset), main="Expression levels of RELN according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="RELN",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="WNT10B",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of WNT10B according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="WNT10B",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="FZD8",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of FZD8 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="FZD8",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="TGFB3",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of TGFB3 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="TGFB3",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="JAG2",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of JAG2 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="JAG2",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="NOTCH3",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of NOTCH3 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="NOTCH3",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="GLI1",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of GLI1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="GLI1",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="TENM4",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of TENM4 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="TENM4",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="BMP1",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of BMP1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="BMP1",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="FZD10",]))~eset$WNTACTIVE, data=exprs(eset), main="Expression levels of FZD10 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="FZD10",]))~eset$WNTACTIVE,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)

boxplot(as.vector(exprs(eset[rownames(eset)=="TERF1",]))~eset$WNT, data=exprs(eset), main="Expression levels of TERF1 according to mutation type", xlab="Mutation type", ylab="Log2 Expression")
stripchart(as.vector(exprs(eset[rownames(eset)=="TERF1",]))~eset$WNT,data=exprs(eset),vertical=T,method="jitter",pch=21,col="black",bg="yellow",add=T)




### ANOVA
fit <- aov(as.vector(exprs(eset[rownames(eset)=="KCNJ5",]))~eset$WNT, data=data.frame(exprs(eset)))
plot(fit)
summary(fit) # display Type I ANOVA table
drop1(fit,~.,test="F") # type III SS and F Tests
TukeyHSD(fit) # where fit comes from aov()
$`eset$mRNA_K5`
> summary(fit) # display Type I ANOVA table
Df Sum Sq Mean Sq F value  Pr(>F)    
eset$WNT     3  77.71  25.902      16 3.9e-08 ***
  Residuals   75 121.42   1.619                    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> drop1(fit,~.,test="F") # type III SS and F Tests
Single term deletions

Model:
  as.vector(exprs(eset[rownames(eset) == "AXIN2", ])) ~ eset$WNT
Df Sum of Sq    RSS    AIC F value    Pr(>F)    
<none>                121.42 41.954                      
eset$WNT  3    77.705 199.12 75.034  15.999 3.905e-08 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> TukeyHSD(fit) # where fit comes from aov()
Tukey multiple comparisons of means
95% family-wise confidence level

Fit: aov(formula = as.vector(exprs(eset[rownames(eset) == "KCNJ5", ])) ~ eset$WNT, data = data.frame(exprs(eset)))

$`eset$WNT`
diff        lwr        upr     p adj
CTNNB1-APC    3.3717720  0.8183301  5.9252139 0.0047147
WT-APC        0.6046828 -1.8052528  3.0146184 0.9119467
ZNRF3-APC     1.3673068 -1.1599453  3.8945588 0.4900293
WT-CTNNB1    -2.7670892 -3.8397493 -1.6944290 0.0000000
ZNRF3-CTNNB1 -2.0044652 -3.3196889 -0.6892416 0.0008200
ZNRF3-WT      0.7626240 -0.2461065  1.7713545 0.2022529

# Class Comparison -- limma
# designmatrix 1: TYPE (CTNNB1 x WT x ZNFR3 x APC)
library(limma)
design1 <- model.matrix(~ 0+pData(eset)$WNTACTIVE)
colnames(design1) <- c("WNT_high","WNT_int","WNT_low")

fit <- lmFit(eset, design1)

contrast.matrix <- makeContrasts(WNT_high-WNT_int, WNT_high-WNT_low, WNT_int-WNT_low, levels=design1)


fit <- eBayes(contrasts.fit(fit, contrast.matrix))

dt <- decideTests(fit, adjust.method="BH", p.value=0.05, lfc=log2(1.3))

vennDiagram(dt)
results <- vennCounts(dt)
vennDiagram(dt,include=c("up","down"),counts.col=c("red","green"),show.include=F)



# Writing tables
## by type
### contrast 1
# Writing tables - corrected by array weight
## by type
### contrast 1 - PS
write.table(topTable(fit, coef=1, adjust="BH", number=Inf, p.value=0.05, lfc=log2(1.3)), file="WNT_high-WNT_int.xls", row.names=TRUE, col.names=TRUE, sep="\t")
write.table(topTable(fit, coef=2, adjust="BH", number=Inf, p.value=0.05, lfc=log2(1.3)), file="WNT_high-WNT_low.xls", row.names=TRUE, col.names=TRUE, sep="\t")
write.table(topTable(fit, coef=3, adjust="BH", number=Inf, p.value=0.05, lfc=log2(1.3)), file="WNT_int-WNT_low.xls", row.names=TRUE, col.names=TRUE, sep="\t")


pData(eset)<-cbind(pData(eset),rownames(pData(eset)))
colnames(pData(eset))[24]<-"Array"
library(maanova)
ACC<-read.madata(datafile=exprs(eset), designfile=pData(eset), header=T, matchDataToDesign=T)


C=matrix(c(0,1,-1,0,0,0,-1,1,1,0,-1,0,0,1,0,-1), ncol=4, byrow=T)
fit.full.mix<-fitmaanova(ACC,formula=~0+WNT,random=~1,covariate=~1)
resiplot(ACC,fit.full.mix)
ftest.all = matest(ACC, fit.full.mix, test.method=c(1,1),shuffle.method="sample", term="WNT", n.perm= 100)
ftest.pair=matest(ACC,fit.full.mix,Contrast=C,term="WNT",n.perm=100)
ftest.all = adjPval(ftest.all, method = "stepup")
ftest.pair = adjPval(ftest.pair, method = "stepup")
summarytable(ftest.pair,method=c("adjPvalperm","Fold.change"), test="F1",outfile="WNT.xls")
volcano(ftest.pair)
idx.fix = volcano(ftest.pair)
cluster.kmean <- macluster(fit.full.mix, term="WNT",idx.gene=idx.fix$comparison1$idx.Fs,what="gene", method="kmean",kmean.ngroups=5, n.perm=100)
con.kmean <- consensus(cluster.kmean, 0.7)
con.kmean$groupname
cluster.hc <- macluster(fit.full.mix, term="WNT",idx.gene=idx.fix$comparison1$idx.Fs,what="sample", method="hc", n.perm=100)
con.hc <- consensus(cluster.hc)


### WNTACTIVE
C=matrix(c(1,-1,0,1,0,-1,0,1,-1), ncol=3, byrow=T)
fit.full.mix<-fitmaanova(ACC,formula=~0+WNTACTIVE,random=~1,covariate=~1)
resiplot(ACC,fit.full.mix)
ftest.all = matest(ACC, fit.full.mix, test.method=c(1,1),shuffle.method="sample", term="WNTACTIVE", n.perm= 100)
ftest.pair=matest(ACC,fit.full.mix,Contrast=C,term="WNTACTIVE",n.perm=100)
ftest.all = adjPval(ftest.all, method = "stepup")
ftest.pair = adjPval(ftest.pair, method = "stepup")
summarytable(ftest.pair,method=c("Pvalperm","Fold.change"), test="F1",outfile="WNT_active.xls")
volcano(ftest.pair)
idx.fix = volcano(ftest.pair)
cluster.kmean <- macluster(fit.full.mix, term="WNTACTIVE",idx.gene=idx.fix$comparison1$idx.Fs,what="gene", method="kmean",kmean.ngroups=6, n.perm=100)
con.kmean <- consensus(cluster.kmean, 0.7)
con.kmean$groupname
cluster.hc <- macluster(fit.full.mix, term="WNTACTIVE",idx.gene=idx.fix$comparison1$idx.Fs,what="sample", method="hc", n.perm=100)
con.hc <- consensus(cluster.hc)

### mRNA_K5
C=matrix(c(1,-1,0,0,0,1,0,-1,0,0,1,0,0,-1,0,1,0,0,0,-1,0,1,-1,0,0,0,1,0,-1,0,0,1,0,0,-1,0,0,1,-1,0,0,0,1,0,-1,0,0,0,1,-1), ncol=5, byrow=T)
fit.full.mix<-fitmaanova(ACC,formula=~0+mRNA_K5,random=~1,covariate=~1)
resiplot(ACC,fit.full.mix)
ftest.all = matest(ACC, fit.full.mix, test.method=c(1,1),shuffle.method="sample", term="mRNA_K5", n.perm= 100)
ftest.pair=matest(ACC,fit.full.mix,Contrast=C,term="mRNA_K5",n.perm=100)
ftest.all = adjPval(ftest.all, method = "stepup")
ftest.pair = adjPval(ftest.pair, method = "stepup")
summarytable(ftest.pair,method=c("Pvalperm","Fold.change"), test="F1",outfile="mRNA_K5.xls")
volcano(ftest.pair)
idx.fix = volcano(ftest.pair)
cluster.kmean <- macluster(fit.full.mix, term="mRNA_K5",idx.gene=idx.fix$comparison1$idx.Fs,what="gene", method="kmean",kmean.ngroups=5, n.perm=100)
con.kmean <- consensus(cluster.kmean, 0.7)
con.kmean$groupname
cluster.hc <- macluster(fit.full.mix, term="mRNA_K5",idx.gene=idx.fix$comparison1$idx.Fs,what="sample", method="hc", n.perm=100)
con.hc <- consensus(cluster.hc)


library(calibrate)

xprs <- data.frame(t(exprs(eset)))

qplot(C9orf84, CSGALNACT1, data=xprs, color=eset$mRNA_K5,size=SHH)

MBNL2
      
qplot(LHX2, AXIN2, data=xprs, color=eset$SCNA_CLS, size=MKI67)

###boxplot
qplot(eset$WNT,IGF2, data=xprs, size=I(1.5), color=eset$WNT, geom=c("boxplot", "jitter"))




textxy(exprs(eset[rownames(eset)=="C9orf84"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="C9orf84"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="C9orf84"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="PCP4"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="PCP4"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="ISM1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="ISM1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="RELN"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="RELN"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="VSNL1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="VSNL1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="CHL1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="CHL1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="LGR5"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="LGR5"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), eset$WNT)

plot(exprs(eset[rownames(eset)=="NKD1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="NKD1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="NOV"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="NOV"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="RHBG"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="RHBG"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="KCNJ5"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="KCNJ5"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="SLC30A10"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="SLC30A10"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="COL26A1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="COL26A1"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="ABCB4"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="ABCB4"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="PDE2A"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="PDE2A"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="MYH7B"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="MYH7B"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="LAMC3"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="LAMC3"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="AFF3"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="AFF3"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="GOLGA7B"]), exprs(eset[rownames(eset)=="CSGALNACT1"]))
textxy(exprs(eset[rownames(eset)=="GOLGA7B"]), exprs(eset[rownames(eset)=="CSGALNACT1"]), colnames(eset))

plot(exprs(eset[rownames(eset)=="GLI1"])[eset$WNT=="WT"], exprs(eset[rownames(eset)=="AXIN2"])[eset$WNT=="WT"])
textxy(exprs(eset[rownames(eset)=="GLI1"])[eset$WNT=="WT"], exprs(eset[rownames(eset)=="AXIN2"])[eset$WNT=="WT"], eset$mRNA_K2[eset$WNT=="WT"])

plot(exprs(eset[rownames(eset)=="CYP11"]), exprs(eset[rownames(eset)=="AXIN2"]))
textxy(exprs(eset[rownames(eset)=="GLI1"]), exprs(eset[rownames(eset)=="AXIN2"]), eset$WNT)
plot(exprs(eset[rownames(eset)=="SHH"]), exprs(eset[rownames(eset)=="AXIN2"]))
textxy(exprs(eset[rownames(eset)=="SHH"]), exprs(eset[rownames(eset)=="AXIN2"]), eset$WNT)

library(STRINGdb)
string_db <- STRINGdb$new(version="9_1", species=9606, score_threshold=0, input_directory="")

gene<-"MBNL2"

#cutoff <- median(exprs(eset[rownames(eset)==gene])[eset$WNT=="CTNNB1"])
#colnames(eset)[eset$WNT=="WT" & exprs(eset)[gene,] > cutoff]

#cutoff<-exprs(eset[rownames(eset)==gene])[eset$WNT=="CTNNB1"][order(exprs(eset[rownames(eset)==gene])[eset$WNT=="CTNNB1"])[1]]
### correlations coefficients


#funccor<-function(x) {cor(exprs(eset)[gene,][eset$WNT=="ZNRF3"],x, method = "pearson")}
funccor<-function(x) {cor(exprs(eset)[gene,],x, method = "pearson")}
#funccortest<-function(x) {cor.test(exprs(eset)[gene,][eset$WNT=="ZNRF3"],x, method = "pearson",conf.level = 0.95,alternative="two.sided",na.rm=TRUE,adjust="BH",alpha=.05)}
funccortest<-function(x) {cor.test(exprs(eset)[gene,],x, method = "pearson",conf.level = 0.95,alternative="two.sided",na.rm=TRUE,adjust="BH",alpha=.05)}
#funccortest<-function(x) {bayes.cor.test(exprs(eset)[gene,][eset$WNT=="ZNRF3"],x, method = "pearson",conf.level = 0.95,alternative="two.sided",na.rm=TRUE,adjust="BH",alpha=.05)}
#cor.gene<-apply(exprs(eset)[,eset$WNT=="ZNRF3"],1,funccor)
cor.gene<-apply(exprs(eset),1,funccor)
#selected<-head(order(cor.gene, decreasing=TRUE),100)
#cor.genetest<-apply(exprs(eset)[,eset$WNT=="ZNRF3"],1,funccortest)
cor.genetest<-apply(exprs(eset),1,funccortest)

pfromlist <- function(x) { unlist(x)[3] }
p<-unlist(lapply(cor.genetest,pfromlist))
corr <- cbind(as.numeric(cor.gene),as.numeric(p))
colnames(corr)<-c("cor.gene","p")
rownames(corr)<-rownames(exprs(eset))
selected <- corr[,2]<0.01 & corr[,1]>0
corr<-corr[selected,]
ordered <- order(corr[,2])
corr <- corr[ordered,]
corrnames<-rownames(corr)
z<-data.frame(corr,corrnames)
colnames(z)<-c("logFC","pvalue","gene")


mapped <- string_db$map(z, "gene", removeUnmappedRows = TRUE)
hits <- mapped$STRING_id[1:200]
z_color <- string_db$add_diff_exp_color(z,logFcColStr="logFC")
payload_id <- string_db$post_payload(mapped$STRING_id,colors=z_color$color)
string_db$plot_network(hits) #,payload_id=payload_id)

string_db$plot_ppi_enrichment(mapped$STRING_id, quiet=TRUE)

### perform gene set enrichment analysis


#funccor<-function(x) {cor(exprs(eset)[gene,],x)}
#cor.gene<-apply(exprs(eset),1,funccor)
#head(sort(cor.gene, decreasing=FALSE),50)

write.table(corr,"NANOG_all.xls", sep="\t")

### building a heatmeap with wnt-specific genes

wntgenes<-exprs(eset)[c("ISM1","PCP4","KCNJ5","PDE2A","LGR5","CDH2","WNT4","NKD1","TCF7","AXIN2","LEF1","BMP4","TERF1","RNF43","ZNRF3"),]

wntgenes<-exprs(eset)[c("NANOG","TERF1","LGR5"),]

mscgenes<-exprs(eset)[c("ANPEP", "ITGB1", "CD44", "ITGA5", "ICAM1", "TFRC", "NT5E", "THY1", "ENG", "VCAM1", "ALCAM"),]


wntgenes3<-exprs(eset)[c("AXIN2","LEF1","TCF7","ISM1","LGR5","TERF1","NANOG","MKI67","TOP2A","TERT"),]


wntgenes2<-wntgenes[,order(comp1$contrib)]

wntgenes3<-wntgenes3[,order(comp1$contrib)]

#wntgenes2<-wntgenes[,eset$WNT=="WT"]

hc2 <- reorder(as.dendrogram(hclust(dist(t(mscgenes),method="euclidean")),method="average"),79:1)

                   
library(gplots)
test <- heatmap.2(mscgenes, col=bluered(75),dendrogram="col",Rowv=F,Colv=hc2,trace="none",scale="row",labCol=pData(eset)[,"WNT"][order(comp1$contrib)],key=F,colsep=c(16,55), sepcolor="green",
                  ColSideColors=c("yellow",rep("green",5),"magenta",rep("green",20),"yellow",rep("green",4),"yellow",rep("green",4),"yellow",rep("green",4),
                                  rep("yellow",2),"green","yellow",rep("green",5),"yellow",rep("green",2),"yellow","green",rep("yellow",2),
                                  "green","yellow","green",rep("yellow",2),rep("green",2),rep("red",9),"magenta","red","green",rep("red",2)),
                  margins=c(5,5),sepwidth=(c(0.25,0.25)))
par(lend=1)
legend("bottomleft",      # location of the legend on the heatmap plot
       legend = c("CTNNB1", "APC", "ZNRF3","No mutation"), # category labels
       col = c("red", "magenta", "yellow","green"),  # color key
       lty= 5,             # line style
       lwd = 8,            # line width
       text.font=1,
       cex=.84
)


test2 <- heatmap.2(wntgenes3, col=bluered(75),dendrogram="none",Rowv=F,Colv=F,trace="none",scale="row",labCol=pData(eset)[,"WNT"][order(comp1$contrib)],key=F, sepcolor="green",
                 ColSideColors=c("yellow",rep("green",5),"magenta",rep("green",20),"yellow",rep("green",4),"yellow",rep("green",4),"yellow",rep("green",4),
                                 rep("yellow",2),"green","yellow",rep("green",5),"yellow",rep("green",2),"yellow","green",rep("yellow",2),
                                 "green","yellow","green",rep("yellow",2),rep("green",2),rep("red",9),"magenta","red","green",rep("red",2)),
                 margins=c(5,5),sepwidth=(c(0.25,0.25)))

par(lend=1)
legend("bottomleft",      # location of the legend on the heatmap plot
       legend = c("CTNNB1", "APC", "ZNRF3","No mutation"), # category labels
       col = c("red", "magenta", "yellow","green"),  # color key
       lty= 5,             # line style
       lwd = 8,            # line width
       text.font=1,
       cex=.84
)



test <- heatmap.2(wntgenes2, col=bluered(75),trace="none",scale="row",labCol=pData(eset)[,"WNT"][order(comp1$contrib)],key=F,colsep=c(27,65), sepcolor="green",sepwidth=(c(0.15,0.15)))
order(comp1$contrib, decreasing=T)

### creating a heatmap using the pachage heatplus
# 1. regular heatmaps
reg1 <- regHeatmap(wntgenes2,col=bluered,dendrogram=list(clustfun=hc))


plot(reg1)

corrdist=function(x) as.dist(1-cor(t(x)))
hclust.avl=function(x) reorder(as.dendrogram(hclust(x,method="average"),79:1))

clust=function(x) {heatmap.2(wntgenes2, col=bluered(75),dendrogram="none",Rowv=F,Colv=F,trace="none",scale="row",labCol=pData(eset)[,"WNT"][order(comp1$contrib)])}

reg2 <- annHeatmap(wntgenes2,
                   dendrogram=list("none",legend=3, col=bluered,annotation=list(data=pData(eset)[,"WNT"]),)
plot(reg2)

reg4 <- annHeatmap2(exprs(eset)[c("ISM1","LGR5","WNT4","NKD1","TCF7","CTNNA2","PDE2A","LEF1","ZNRF3","RNF43","WIF1","BMP4"),], legend=3, Col=list(status="hide"),Row=list(status="hide")),col=heat.colors,annotation=list(pData(eset)[,"WNT"]))
plot(reg4)

reg5 <- annHeatmap2(exprs(eset)[c("ISM1","LGR5","WNT4","NKD1","TCF7","CTNNA2","PDE2A","LEF1","ZNRF3","RNF43","WIF1","BMP4"),], dendrogram=list(heatmap.2(data="none"), annotation=list(data=pData(eset)[,"WNT"]))
plot(reg5)

library(reshape2)
library(ggplot2)
library(RColorBrewer)

######

wntframe <- data.frame(wntgenes)

head(wntframe)









myPalette <- colorRampPalette(brewer.pal(9, "Spectral"))
wntmelt<-melt(wntgenes)
zp1 <- ggplot(wntmelt, aes(x = X2, y = X1, fill = value))
zp1 <- zp1 + geom_tile()
p1 <- zp1 + scale_fill_gradientn(colours = myPalette(100))
zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
zp1 <- zp1 + coord_equal()
zp1 <- zp1 + theme_bw()
print(zp1)
library(lattice)
heatmap(wntgenes, Rowv = NA, Colv = NA, scale = "none", col = myPalette)
