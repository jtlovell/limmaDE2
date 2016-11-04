## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=6,echo=TRUE, warning=FALSE, message=FALSE)

## ----install limmaDE2, eval = F------------------------------------------
#  library(devtools)
#  install_github("jtlovell/limmaDE2")

## ----load packages-------------------------------------------------------
library("limmaDE2")
library("ggplot2")
library("reshape2")
library("SimSeq")
library("DESeq2")

## ----load data-----------------------------------------------------------
data(kidney) # from simseq
counts<-kidney$counts
counts<-counts[sample(1:nrow(counts),1000),]
info<-with(kidney, 
           data.frame(id = paste(replic, treatment, sep = "_"), 
                      rep=replic, 
                      Treatment=ifelse(treatment == "Tumor","tumor","cntr"), 
                      stringsAsFactors=F))
colnames(counts)<-info$id

## ----organize data-------------------------------------------------------
group.a<-c("4619", "4712", "4863", "4865", "5452", "5453", "5454", "5455",
          "5456", "5457","5458", "5461", "5462", "5463", "5465", "5467",
          "5468", "5469", "5470", "5549","5552", "5580", "5641", "5672",
          "5689", "5701", "5703", "5706", "5989", "6088")
info$group<-ifelse(info$rep %in% group.a, "a","b")
with(info, table(group, Treatment))

info$trt.grp<-with(info, paste(Treatment, group, sep="_"))


## ------------------------------------------------------------------------
head(info)

## ------------------------------------------------------------------------
counts[1:10,1:3]

## ------------------------------------------------------------------------
info$Treatment <- factor(info$Treatment,
                           levels = c("cntr", "tumor"))
info$group <- factor(info$group,
                           levels = c("a", "b"))

## ----basic pipe----------------------------------------------------------
stats <- pipeLIMMA(counts = counts, 
                   info = info, 
                   block = NULL, 
                   formula = "~ Treatment")
lmStats<-stats$stats
voom<-stats$voom$E

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(stats$stats[1:6,1:6])

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(stats$voom[["E"]][1:6,1:6])

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(stats$fstats[1:6,])

## ----pipe w block--------------------------------------------------------
info$block <- rep(1:2,each=nrow(info)/2)
stats.block <- pipeLIMMA(counts = counts, 
                   info = info, 
                   block = info$block,
                   formula = "~ Treatment")

## ----pipe factorial------------------------------------------------------
stats.factorial <- pipeLIMMA(counts = counts, 
                   info = info, 
                   block = NULL, 
                   formula = "~ Treatment + group + Treatment*group")

## ------------------------------------------------------------------------
design <- with(info, model.matrix(~ 0 + trt.grp))
colnames(design)<-gsub("trt.grp","",colnames(design))
head(design)

## ------------------------------------------------------------------------
contrast.matrix <- makeContrasts(
  tumor_a - cntr_a , 
  tumor_b - cntr_b,
  levels = design)
head(contrast.matrix)

## ----contrast pipe-------------------------------------------------------
stats <- pipeLIMMA(counts = counts, 
                   info = info, 
                   block = NULL, 
                   design = design, 
                   contrast.matrix = contrast.matrix)
stats.contrasts <- stats$stats

## ------------------------------------------------------------------------
sigs <- makeBinarySig(stats.contrasts, what = "Qvalue", alpha = 0.05)

## ------------------------------------------------------------------------
counts2Venn(x=sigs, cols=c(1,2), names=c("in.grpA","in.grpB"),
   colors=c("blue","red"),type="limma", legx=-3.3,legy=-3)

## ------------------------------------------------------------------------
counts2Venn(x=sigs, cols=c(1,2), names=c("in.grpA","in.grpB"),
   colors=c("blue","red"),type="Euler", legx=-3.3,legy=-3)

## ------------------------------------------------------------------------
with(lmStats, volcanoPlot(pval = ebayesPvalue_Treatmenttumor,
                          lfc = Treatmenttumor_logFC,
                          sig = ebayesQvalue_Treatmenttumor,
                          main = "no tumor vs. tumor Volcano Plot", 
                          xlab = "tumor - no tumor Log2 Fold Change",
                          bty = "n", legpos = "top", leginset = c(0,-.1)))

## ------------------------------------------------------------------------
sigs <- data.frame(makeBinarySig(stats.contrasts, what = "Qvalue", alpha = 0.05))
names(sigs)<-c("sig.a","sig.b")
sigs$sign.a<-sign(stats.contrasts$tumor_a...cntr_a_logFC)
sigs$sign.b<-sign(stats.contrasts$tumor_b...cntr_b_logFC)

cols<-with(sigs, ifelse(sig.a + sig.b == 0,  
                        "grey",
                        ifelse(sig.a + sig.b == 2 & sign.a*sign.b == 1, 
                               "pink",
                               ifelse(sig.a + sig.b == 2 & sign.a*sign.b == -1,
                                      "cornflowerblue",
                                      ifelse(sig.a == 1, "darkblue", "darkred")))))

## ----pac-----------------------------------------------------------------
with(stats.contrasts, volcanoPair(lfc1 = tumor_a...cntr_a_logFC,
                                  lfc2 = tumor_b...cntr_b_logFC,
                                  pt.col = cols, pt.pch = 19, pt.cex=.5,
                                  xlab = "Tumor - control LFC (group A)",
                                  ylab = "Tumor - control LFC (group B)"))

## ------------------------------------------------------------------------
pca12 <- counts2PCA(counts=voom, info = info, ids = info$id, pcas2return = 6)
pca12.var <- pca12[[2]]
pca12 <- pca12[[1]]
gcols <- c("darkblue", "blue", "gold", "green", "pink", "red")
ggplot(pca12, aes(x = PC1, y = PC2, col = group, shape = Treatment)) +
  geom_vline(xintercept = 0, lty = 2)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  labs(x = paste("PCA1 (", pca12.var[1],"%)", sep = ""),
       y = paste("PCA2 (", pca12.var[2],"%)", sep = ""),
       title = "PCA analysis of voom-normalized expression")

