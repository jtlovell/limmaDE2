# change to your local directory
setwd("/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/run_full_analysis")

#clear the environment
rm(list=ls())

#install the package with R functions
library(devtools)
install_github("jtlovell/limmaDE2")
library(limmaDE2)
# Note 1. to load this package, you need to have the following packages installed:
# limma, SimSeq
# Install by:
# source("https://bioconductor.org/biocLite.R")
# biocLite() # only do this if you have not installed the base bioconductor packages -
#     ### can take a while
# biocLite(c("limma","SimSeq"))

# Note 2. To use some of the functions in this package you need a few other packages installed:
# ggplot2, edgeR, venneuler, qvalue, qdap, reshape2, dendroextras, Heatplus, gdata
# install by:
# biocLite(c("Heatplus", "qvalue","edgeR"))
# install.packages("ggplot2", "venneuler", "qdap", "reshape2", "dendroextras", "gdata")

#install the package with data
install_github("jtlovell/PhyGenomicsData.PV2016")
library(PhyGenomicsData.PV2016)

# Run just the 2012 data analysis:
### for conviencence, rename info and counts
info<-info12
counts<-counts12

#################################
# Part 1: Full model with Treatment
#################################
## Run LIMMA, with subplot treated as repeated measures of the same plant.
### The only response is treatment. Sub_Block is the sub-block in the split plot design, containing two AP13
###    plants. The duplicate correlation and blocking in LIMMA should account for correlations between these plants
### We do not estimate intercept, so that the model only looks at the effects of Treatment, and the F-test reflects this
## Run the modeling pipeline
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ Treatment")
v<-stats$voom[["E"]]
stats.model1<-stats$stats

#################################
# Part 2: All pairwise contrasts
#################################
## Generate a design matrix that specifies all relavant contrasts
f<-info$Treatment
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(
  f25th-flow, fmean-flow, fambient-flow, f75th-flow, fhigh-flow, fmean-f25th, fambient-f25th, fhigh-f25th,
  fambient-fmean, f75th-fmean, fhigh-fmean, f75th-fambient, fhigh-fambient, fhigh-f75th, levels=design)
stats<-pipeLIMMA(counts=counts, info=info, block=NULL, design=design, contrast.matrix=contrast.matrix)
stats.model1.contrasts<-stats$stats

## Look at some of the effects, for example, side-by-side volcano plots
par.original <- par(no.readonly=TRUE)
par(mfrow=c(2,2))
with(stats.model1.contrasts,
     volcanoPlot(pval=ebayesPvalue_fhigh...flow, lfc=fhigh...flow_logFC,
                 sig=ebayesQvalue_fhigh...flow, main="high vs. low",
                 ylim=c(0,7), xlim=c(-6.1,6.1), legpos=c(-.6,0)))
with(stats.model1.contrasts,
     volcanoPlot(pval=ebayesPvalue_f75th...flow,  lfc=f75th...flow_logFC,
                 sig=ebayesQvalue_f75th...flow, main="75th vs. low",
                 ylim=c(0,7), xlim=c(-6.1,6.1), legpos=c(-.6,0)))
with(stats.model1.contrasts,
     volcanoPlot(pval=ebayesPvalue_fmean...flow, lfc=fmean...flow_logFC,
                 sig=ebayesQvalue_fmean...flow, main="mean vs. low",
                 ylim=c(0,7), xlim=c(-6.1,6.1), legpos=c(-.6,0)))
with(stats.model1.contrasts,
     volcanoPlot(pval=ebayesPvalue_f25th...flow, lfc=f25th...flow_logFC,
                 sig=ebayesQvalue_f25th...flow, main="25th vs. low",
                 ylim=c(0,7), xlim=c(-6.1,6.1), legpos=c(-.6,0)))
par(par.original)

## Also, plot the effects across treatments for the top 20 genes
tp.trt<-topGeneRxn(v=v, info=info,
                     sig=stats.model1.contrasts$Fqvalue,
                     pointcols=NULL, paletteChoice="RdBu",
                     xdat="Treatment")
tp.mdwp<-topGeneRxn(v=v, info=info,
                     sig=stats.model1.contrasts$Fqvalue,
                     pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                     xdat="MDWP", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)")


#################################
# Part 3: Run PCA
#################################
pca<-voom2PCA(v=v, info=info, ids=info$ID)
ggplot(pca, aes(x=PC1, y=PC2, col=Treatment, alpha=Treatment))+
  theme_bw()+geom_point(size=4)

#################################
# Part 4: Subset to lines with physiology data
#################################
phys.lines<-info$ID[!is.na(info$order)]
counts<-counts[,phys.lines]
info<-info[info$ID %in% phys.lines,]
info$Treatment<-factor(info$Treatment, levels=c("low","mean","high"))

#################################
# Part 5: Run model with MDWP as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ MDWP")
v<-stats$voom[["E"]]
stats.allests.mdwp<-stats$stats
par(mfrow=c(1,1))
with(stats.allests.mdwp,
     volcanoPlot(pval=ebayesPvalue_MDWP, lfc=MDWP_logFC,
                 sig=ebayesQvalue_MDWP, main="high vs. low"))
# plot the top genes for mdwp
tp.mdwp<-topGeneRxn(v=v, info=info,
                    sig=stats.allests.mdwp$ebayesQvalue_MDWP,
                    pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                    xdat="MDWP", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)")
# pairwise plot of treatment (hvl and mdwp)
volcanoPair(lfc1=stats.allests.mdwp$MDWP_logFC,
            lfc2=stats.model1.contrasts$fhigh...flow_logFC,
            sig1=stats.allests.mdwp$ebayesQvalue_MDWP,
            sig2=stats.model1.contrasts$ebayesQvalue_fhigh...flow,
            xlab="mdwp Log2 Fold Change", ylab="treatment high vs. low log2 fold change")
# pull out strongest genes that have mdwp effect but not treatment
wp.genes<-stats.allests.mdwp$gene[stats.allests.mdwp$ebayesQvalue_MDWP<=0.01 & stats.model1.contrasts$ebayesQvalue_fhigh...flow>=0.05]
tp.mdwp<-topGeneRxn(v=v, info=info,
                    geneIDs=wp.genes,
                    sig=stats.allests.mdwp$ebayesQvalue_MDWP,
                    pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                    xdat="MDWP", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)",
                    title="effect of genes with strong MDWP effects, \n but weak treatment effects")
# make venn diagram of the number of gene effected by MDWP and treatment
mdwp.trt.Qvals<-data.frame(MDWP.qvalue=stats.allests.mdwp$ebayesQvalue_MDWP, trt.qvalue=stats.model1.contrasts$ebayesQvalue_fhigh...flow)
sig.q.05<-makeBinarySig(mdwp.trt.Qvals, what="qvalue", alpha=0.05)
counts2Venn(x=sig.q.05, cols=c(1,2), names=c("MDWP","high vs. low"), colors=c("darkblue","red"),
            main="Comparison of models", type="limma", legy=-3.2, legx=-3.5)

#################################
# Part 6: Run model with Sampling Order as the predictor
#################################
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ order")
stats.allests.order<-stats$stats
with(stats.allests.order,
     volcanoPlot(pval=ebayesPvalue_order, lfc=order_logFC,
                 sig=ebayesQvalue_order, main="high vs. low"))
tp.mdwp<-topGeneRxn(v=v, info=info,
                    sig=stats.allests.order$ebayesQvalue_order,
                    pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                    xdat="order", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)")

# save(stats.fullmodel, stats.allests, lim.contrasts, pca, v,
#      stats.fullmodel.mdwp, stats.allests.mdwp, stats.fullmodel.order, stats.allests.order,stats.fullmodel.subtrt, stats.allests.subtrt,
#      file="/Users/John/Desktop/dropbox/Switchgrass_PlantPhys/stats_output/tempe2012_allstats.RData")

