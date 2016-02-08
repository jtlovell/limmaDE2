# change to your local directory
setwd("/Users/JLovell/dropbox/Switchgrass_PlantPhys/run_full_analysis/y2012")

#clear the environment
rm(list=ls())

#install the package with R functions
library(devtools)
install_github("jtlovell/limmaDE2")
library(limmaDE2)
library(ggplot2)
library(RColorBrewer)
library(Heatplus)
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
data(info12)
data(counts12)
info<-info12
counts<-counts12

#################################
# Part 1: Full model with Treatment
#################################
## 1.1: Summary of experimental design
pdf("experimentalDesign_plots.pdf", height=4, width=6)
png("experimentalDesign_plots.png", width = 300, height = 240)
ggplot(info[!is.na(info$PDWP),], aes(x=PDWP, y=MDWP, col=Treatment))+geom_point(size=4)+ theme_bw() +
  geom_abline(slope=1, intercept=0, lty=2)+
  scale_color_manual(values=c("darkred","forestgreen","skyblue"))+
  scale_x_continuous("Pre-dawn Leaf Water Potential (MPa)")+
  scale_y_continuous("Mid Day Leaf Water Potential (MPa)")+
  theme(panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  ggtitle("Effect of sampling order \non measured leaf water potential")
dev.off()

################
## 1.2: Run the limma pipeline
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ Treatment")
v<-stats$voom[["E"]]
stats.model1<-stats$stats

################
## 1.3: Heatmap of highly significant genes
pdf("maineffect_heatmaps.pdf", height=4, width=6)
mean.v<-means4heatmap(v[stats.model1$Fqvalue<=0.01,], grps=info$Treatment)
v2<-data.matrix(mean.v); rownames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red'))))
v2<-data.matrix(v[stats.model1$Fqvalue<=0.01,]); rownames(v2)<-NULL; colnames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red')),
                 ann = list(Col = list(data = info[,c("Treatment")]))))
dev.off()

png("maineffect_mean_heatmap.png", width = 480, height = 400)
mean.v<-means4heatmap(v[stats.model1$Fqvalue<=0.01,], grps=info$Treatment)
v2<-data.matrix(mean.v); rownames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red'))))
dev.off()
png("maineffect_ind_heatmap.png", width = 480, height = 400)
v2<-data.matrix(v[stats.model1$Fqvalue<=0.01,]); rownames(v2)<-NULL; colnames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red')),
                 ann = list(Col = list(data = info[,c("Treatment")]))))
dev.off()

#################################
# Part 2: All pairwise contrasts
#################################
## 2.1. Make contrast matrix
f<-info$Treatment
design <- model.matrix(~0+f)
contrast.matrix<-makeContrasts(
  f25th-flow, fmean-flow, fambient-flow, f75th-flow, fhigh-flow, fmean-f25th, fambient-f25th, fhigh-f25th,
  fambient-fmean, f75th-fmean, fhigh-fmean, f75th-fambient, fhigh-fambient, fhigh-f75th, levels=design)

################
## 2.2: Run the limma pipeline
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, design=design, contrast.matrix=contrast.matrix)
stats.model1.contrasts<-stats$stats

################
## 2.3 Look at some of the effects, for example, side-by-side volcano plots
pdf("treatmentContrast_volcanoPlots.pdf", height=7, width=10)
png("treatmentContrast_volcanoPlots.png", width = 480, height = 400)
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
dev.off()

################
## 2.4 Also, plot the effects across treatments for the top 20 genes
pdf("rxnNorms_topContrastGenes.pdf")
png("rxnNorms_topContrastGenes_boxplots.png")
tp.trt<-topGeneRxn(v=v, info=info,
                     sig=stats.model1.contrasts$Fqvalue,
                     pointcols=NULL, paletteChoice="RdBu",
                     xdat="Treatment")
png("rxnNorms_topContrastGenes_scatterplotsMDWP.png")
tp.mdwp<-topGeneRxn(v=v, info=info,
                     sig=stats.model1.contrasts$Fqvalue,
                     pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                     xdat="MDWP", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)")
dev.off()

#################################
# Part 3: Run PCA
#################################
## 3.1: Run the PCA
pdf("PCA_allgenes_bytreatment.pdf")
png("PCA_allgenes_bytreatment_pairs.png")
pca<-voom2PCA(v=v, info=info, ids=info$ID, pcas2return=6,
              plot.cols=brewer.pal(name="RdBu",n=6)[info$Treatment],
              main="scatterplot matrix of Principal components, \n split by Treatment", pch=19, cex=.8)
dev.off()
################
## 3.2: Plot the PCA
png("PCA_allgenes_bytreatment_ggplot.png", width = 480, height = 400)
ggplot(pca, aes(x=PC1, y=PC2, col=Treatment))+
  stat_ellipse(level=.50)+ggtitle("PCA plot with 50% CI ellipses")+
  geom_hline(yintercept=0, lty=2)+ geom_vline(xintercept=0,lty=2)+
  theme_bw()+geom_point(size=4)+scale_color_manual(values=brewer.pal(n=6,name="RdBu"))+
  theme(panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank())+
  scale_x_continuous("Principal component axis #1 (13.3%)")+
  scale_y_continuous("Principal component axis #2 (11.6%)")
dev.off()

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
################
## 5.1: Experimental design of subset data - physiology
library(reshape2)
pca2<-melt(pca, measure.var=c("PC1","PC2"))
pdf("PCA_allgenes_mdpw_response.pdf", width=6, height=4)
png("PCA_allgenes_mdpw_response.png", width = 480, height = 300)
ggplot(pca2, aes(x=MDWP, y=value))+geom_point(aes(col=Treatment), size=4) + theme_bw() +
  stat_smooth(lty=2, se=F, method="lm", col="black")+
  scale_color_manual(values=brewer.pal(n=6,name="RdBu"))+ facet_grid(.~variable)+
  scale_y_continuous("Principal Component Axis") +
  scale_x_continuous("Mid Day Leaf Water Potential (MPa)") +
  ggtitle("DESeq2 PCA vs. Leaf Water Potential")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
dev.off()
################
## 5.2: Run the limma pipeline
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ MDWP")
v<-stats$voom[["E"]]
stats.mdwp<-stats$stats

################
## 5.3: Make some volcano plots
pdf("mdwp_volcanoPlot.pdf")
png("mdwp_volcanoPlot.png", width = 360, height = 300)
par(mfrow=c(1,1))
with(stats.mdwp,
     volcanoPlot(pval=ebayesPvalue_MDWP, lfc=MDWP_logFC,
                 sig=ebayesQvalue_MDWP, main="high vs. low"))
dev.off()

################
## 5.4: plot the top genes for mdwp
pdf("rxnNorms_topMDWPGenes.pdf")
png("rxnNorms_topMDWPGenes.png")
tp.mdwp<-topGeneRxn(v=v, info=info,
                    sig=stats.mdwp$ebayesQvalue_MDWP,
                    pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                    xdat="MDWP", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)")
dev.off()

################
## 5.5: pairwise plot of treatment (hvl and mdwp)
pdf("pairwiseVolcanoPlot_mdwp.high_v_lowContrast.pdf")
png("pairwiseVolcanoPlot_mdwp.high_v_lowContrast.png", width = 320, height = 300)
volcanoPair(lfc1=stats.mdwp$MDWP_logFC,
            lfc2=stats.model1.contrasts$fhigh...flow_logFC,
            sig1=stats.mdwp$ebayesQvalue_MDWP,
            sig2=stats.model1.contrasts$ebayesQvalue_fhigh...flow,
            xlab="mdwp Log2 Fold Change", ylab="treatment high vs. low log2 fold change")
dev.off()
# pull out strongest genes that have mdwp effect but not treatment

################
## 5.6: rxn norms of mdwp genes without treatment effects
wp.genes<-stats.mdwp$gene[stats.mdwp$ebayesQvalue_MDWP<=0.025 & stats.model1.contrasts$ebayesQvalue_fhigh...flow>=0.1]
pdf("rxnNorms_topMDWPGenes_but_noTrtEffect.pdf")
png("rxnNorms_topMDWPGenes_but_noTrtEffect.png")
tp.mdwp<-topGeneRxn(v=v, info=info,
                    geneIDs=wp.genes,
                    sig=stats.mdwp$ebayesQvalue_MDWP,
                    pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                    xdat="MDWP", coldat="Treatment", xlab="Mid-day Leaf Water Potential (MPa)",
                    title="effect of genes with strong MDWP effects, \n but weak treatment effects")
dev.off()
################
## 5.7: make venn diagram of the number of gene effected by MDWP and treatment
pdf("venn_trtContrast_vs_MDWP.pdf")
png("venn_trtContrast_vs_MDWP.png")
mdwp.trt.Qvals<-data.frame(MDWP.qvalue=stats.mdwp$ebayesQvalue_MDWP, trt.qvalue=stats.model1.contrasts$ebayesQvalue_fhigh...flow)
sig.q.05<-makeBinarySig(mdwp.trt.Qvals, what="qvalue", alpha=0.05)
counts2Venn(x=sig.q.05, cols=c(1,2), names=c("MDWP","high vs. low"), colors=c("darkblue","red"),
            main="Comparison of models", type="limma", legy=3.2, legx=3.5)
dev.off()

################
## 5.7: make heatmaps
pdf("mdwp_effect_heatmaps.pdf")
png("mdwp_effect_heatmaps.png", width = 480, height = 400)
v2<-data.matrix(v[stats.mdwp$ebayesQvalue_MDWP<=0.05,]); rownames(v2)<-NULL; colnames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red')),
                 ann = list(Col = list(data = info[,c("Treatment","MDWP")]))))
dev.off()

#################################
# Part 6: Run model with Sampling Order as the predictor
#################################

## 6.1: Experimental design effect of sampling order
pdf("mdwp_samplingTime_effect.pdf")
png("mdwp_samplingTime_effect.png", width = 480, height = 400)
ggplot(info[!is.na(info$PDWP),], aes(x=order, y=MDWP, col=Treatment))+geom_point(size=4)+ theme_bw() +
  stat_smooth(method="lm", se=F, lty=2, alpha=.5, lwd=.5)+
  scale_color_manual(values=c("darkred","forestgreen","cornflowerblue"))+
  scale_x_continuous("Sampling order (first sampling @ 11:00, last @ 13:00)")+
  scale_y_continuous("Mid Day Leaf Water Potential (MPa)")+
  ggtitle("Effect of sampling order on measured leaf water potential")+
  theme(panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank())
dev.off()
################
## 6.2: Limma pipeline for effect of sampling order
stats<-pipeLIMMA(counts=counts, info=info, block=info$Sub_Block, formula="~ order")
stats.order<-stats$stats

################
## 6.3: Volcano plot of sampling order
# pdf("volcanoPlot_order.pdf")
# with(stats.order,
#      volcanoPlot(pval=ebayesPvalue_order, lfc=order_logFC,
#                  sig=ebayesQvalue_order, main="high vs. low"))
# dev.off()

################
## 6.4: Sampling order effect by top genes
pdf("rxn_norms_top20orderGenes.pdf")
png("rxn_norms_top20orderGenes.png")
tp.mdwp<-topGeneRxn(v=v, info=info,
                    sig=stats.order$ebayesQvalue_order,
                    pointcols=c("darkred","forestgreen","skyblue"), paletteChoice=NULL,
                    xdat="order", coldat="Treatment", xlab="Order of sample collection \n 1st @ 11:00, last @ 13:00")
dev.off()

################
## 6.5: Heatmap
pdf("order_effect_heatmaps.pdf")
png("order_effect_heatmaps_alpha0.05.png", width = 480, height = 400)
v2<-data.matrix(v[stats.order$ebayesQvalue_order<=0.05,]); rownames(v2)<-NULL; colnames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red')),
                 ann = list(Col = list(data = info[,c("Treatment","order")]))))
dev.off()
png("order_effect_heatmaps_alpha0.1.png", width = 480, height = 400)
v2<-data.matrix(v[stats.order$ebayesQvalue_order<=0.1,]); rownames(v2)<-NULL; colnames(v2)<-NULL
plot(annHeatmap2(v2, scale="row",
                 col = colorRampPalette(c('dark blue','white','dark red')),
                 ann = list(Col = list(data = info[,c("Treatment","order")]))))
dev.off()

#################################
# Part 7: Multivariate analysis of expression among any genes that are differentially expresseed
#################################
# make a set of vectors of genes with various DE
genes.treatment<-as.character(stats.model1$gene[stats.model1$Fqvalue<=0.05])
genes.high_v_low<-as.character(stats.model1.contrasts$gene[stats.model1.contrasts$ebayesQvalue_fhigh...flow<=0.05])
genes.mdwp<-as.character(stats.mdwp$gene[stats.mdwp$ebayesQvalue_MDWP<=0.05])
genes.order<-as.character(stats.order$gene[stats.order$ebayesQvalue_order<=0.05])
genes.random10k<-as.character(stats.order$gene[sample(nrow(stats.order),10000)])
genes.allsig<-unique(c(genes.treatment,genes.high_v_low,genes.mdwp,genes.order))

pdf("cca_allgenes_withSigEffect.pdf")
png("cca_allgenes_withSigEffect.png")
cca.allsig<-de2cca(info=info, counts=counts[rownames(counts) %in% genes.allsig,], formula="Treatment + MDWP + order")
dev.off()
pdf("cca_random10kgenes.pdf")
png("cca_random10kgenes.png")
cca.rand<-de2cca(info=info, counts=counts[rownames(counts) %in% genes.random10k,], formula="Treatment + MDWP + order")
dev.off()

#################################
# Part 7: Run GO enrichment analysis for a few sets of significant genes
#################################
data(geneID2GO)
geneNames <- names(geneID2GO)

hvl_go<-pipeTopGo(geneID2GO=geneID2GO, genes.of.interest=genes.high_v_low)
hvl_go.stats<-hvl_go$stats[as.numeric(hvl_go$stats$classicFisher)<1,]
hvl_go.stats$fdr.pval<-p.adjust(hvl_go.stats$classicFisher, method="BH")

order.alpha_0.2<-as.character(stats.order$gene[stats.order$ebayesQvalue_order<=0.2])
order_go<-pipeTopGo(geneID2GO=geneID2GO, genes.of.interest=order.alpha_0.2, toGrep="circadian")
order_go.stats<-order_go$stats[as.numeric(order_go$stats$classicFisher)<1,]
order_go.stats$fdr.pval<-p.adjust(order_go.stats$classicFisher, method="BH")

printGraph(hvl_go$GOdata, hvl_go$resultFisher, firstSigNodes = 1, fn.prefix = "hvl_go", useInfo = "def", pdfSW = TRUE)

save(stats.model1, stats.model1.contrasts, pca, v,
     stats.mdwp, stats.order,
     genes.treatment, genes.high_v_low, genes.mdwp, genes.order, genes.allsig,
     cca.allsig, cca.rand,
     order_go.stats, hvl_go.stats,
     file="temple2012_allstatsV2.RData")

